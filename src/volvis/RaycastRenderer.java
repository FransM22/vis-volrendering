/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package volvis;

import com.jogamp.opengl.GL;
import com.jogamp.opengl.GL2;
import com.jogamp.opengl.util.texture.Texture;
import com.jogamp.opengl.util.texture.awt.AWTTextureIO;
import gui.RaycastRendererPanel;
import gui.TransferFunction2DEditor;
import gui.TransferFunctionEditor;
import java.awt.geom.AffineTransform;
import java.awt.image.AffineTransformOp;
import java.awt.image.BufferedImage;
import java.util.Arrays;
import util.TFChangeListener;
import util.VectorMath;
import volume.GradientVolume;
import volume.Volume;
import volume.VoxelGradient;
//import java.util.Arrays;

/**
 *
 * @author michel
 */
public class RaycastRenderer extends Renderer implements TFChangeListener {

    private Volume volume = null;
    private GradientVolume gradients = null;
    RaycastRendererPanel panel;
    TransferFunction tFunc;
    TransferFunctionEditor tfEditor;
    TransferFunction2DEditor tfEditor2D;
    int rayFunction;
    
    // construct func
    public RaycastRenderer() {
        panel = new RaycastRendererPanel(this);
        panel.setSpeedLabel("0");
    }

    public void setVolume(Volume vol) {
        System.out.println("Assigning volume");
        volume = vol;

        System.out.println("Computing gradients");
        gradients = new GradientVolume(vol);

        // set up image for storing the resulting rendering
        // the image width and height are equal to the length of the volume diagonal
        int imageSize = (int) Math.floor(Math.sqrt(vol.getDimX() * vol.getDimX() + vol.getDimY() * vol.getDimY()
                + vol.getDimZ() * vol.getDimZ()));
        if (imageSize % 2 != 0) {
            imageSize = imageSize + 1;
        }
        image = new BufferedImage(imageSize, imageSize, BufferedImage.TYPE_INT_ARGB);
        
        // An image on which the volume may be rendered. It is scaled to the image field
        createNativeImage();
        
        // create a standard TF where lowest intensity maps to black, the highest to white, and opacity increases
        // linearly from 0.0 to 1.0 over the intensity range
        tFunc = new TransferFunction(volume.getMinimum(), volume.getMaximum());
        
        // uncomment this to initialize the TF with good starting values for the orange dataset 
        //tFunc.setTestFunc();
        
        
        tFunc.addTFChangeListener(this);
        tfEditor = new TransferFunctionEditor(tFunc, volume.getHistogram());
        
        tfEditor2D = new TransferFunction2DEditor(volume, gradients);
        tfEditor2D.addTFChangeListener(this);

        System.out.println("Finished initialization of RaycastRenderer");
    }

    public RaycastRendererPanel getPanel() {
        return panel;
    }

    public TransferFunction2DEditor getTF2DPanel() {
        return tfEditor2D;
    }
    
    public TransferFunctionEditor getTFPanel() {
        return tfEditor;
    }
     

    short getVoxel(double[] coord) {
        if (coord[0] < 0 || coord[0] > volume.getDimX() || coord[1] < 0 || coord[1] > volume.getDimY()
                || coord[2] < 0 || coord[2] > volume.getDimZ()) {
            return 0;
        }
        return (short) Math.floor(volume.getInterpolated(coord[0], coord[1], coord[2]));
    }


    void slicer(double[] viewMatrix) {
        clearImage(nativeImage);

        // vector uVec and vVec define a plane through the origin, 
        // perpendicular to the view vector viewVec
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);

        // image is square
        int imageCenter = nativeImage.getWidth() / 2;

        double[] pixelCoord = new double[3];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);

        // sample on a plane through the origin of the volume data
        double max = volume.getMaximum();
        TFColor voxelColor = new TFColor();

        
        for (int j = 0; j < nativeImage.getHeight(); j++) {
            for (int i = 0; i < nativeImage.getWidth(); i++) {
                pixelCoord[0] = uVec[0] * (i - imageCenter) * (1/renderScale) + vVec[0] * (j - imageCenter) * (1/renderScale)
                        + volumeCenter[0] + viewVec[0] * sampleDepth;
                pixelCoord[1] = uVec[1] * (i - imageCenter) * (1/renderScale) + vVec[1] * (j - imageCenter) * (1/renderScale)
                        + volumeCenter[1] + viewVec[1] * sampleDepth;
                pixelCoord[2] = uVec[2] * (i - imageCenter) * (1/renderScale) + vVec[2] * (j - imageCenter) * (1/renderScale)
                        + volumeCenter[2] + viewVec[2] * sampleDepth;

                int val = getVoxel(pixelCoord);
                
                // Map the intensity to a grey value by linear scaling
                voxelColor.r = val/max;
                voxelColor.g = voxelColor.r;
                voxelColor.b = voxelColor.r;
                voxelColor.a = val > 0 ? 1.0 : 0.0;  // this makes intensity 0 completely transparent and the rest opaque
                // Alternatively, apply the transfer function to obtain a color
                //voxelColor = tFunc.getColor(val);
                
                
                // BufferedImage expects a pixel color packed as ARGB in an int
                int c_alpha = voxelColor.a <= 1.0 ? (int) Math.floor(voxelColor.a * 255) : 255;
                int c_red = voxelColor.r <= 1.0 ? (int) Math.floor(voxelColor.r * 255) : 255;
                int c_green = voxelColor.g <= 1.0 ? (int) Math.floor(voxelColor.g * 255) : 255;
                int c_blue = voxelColor.b <= 1.0 ? (int) Math.floor(voxelColor.b * 255) : 255;
                int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
                nativeImage.setRGB(i, j, pixelColor);
            }
        }
        
        scaleImageTo(nativeImage, image);
    }
    
    // viewVec, coordinate on view plane
    double[] getRayDepth(double[] viewVec, double[] coordOnPlane){
        int[][] border = new int[3][2];
        border[0][0] = 0; // x
        border[0][1] = volume.getDimX();
        border[1][0] = 0; // y
        border[1][1] = volume.getDimY();
        border[2][0] = 0; // z
        border[2][1] = volume.getDimZ();
        
        double fDepth = Double.POSITIVE_INFINITY;
        double bDepth = Double.NEGATIVE_INFINITY;
        
        // traverse over the 3 directions
        for (int d = 0; d < 3; d++){
            if (viewVec[d] != 0){
                // one of dp and dn is negative
                double dp = (border[d][0] - coordOnPlane[d]) / viewVec[d];
                double dn = (border[d][1] - coordOnPlane[d]) / viewVec[d];
                if ((dp > 0 && dn > 0) || (dp < 0 && dn < 0)){
                    // outside the box, no samples
                    return new double[]{0, 0};
                }
                
                if (dp < 0 && dn > 0){
                    double tmp = dp;
                    dp = dn;
                    dn = tmp;
                }
                
                if (dp < fDepth)
                    fDepth = dp;
                
                if (dn > bDepth)
                    bDepth = dn;
            }
            
        }
        double[] fbDepth = new double[] {fDepth, bDepth};
        return fbDepth;
}
    
    void mip(double[] viewMatrix){
        clearImage(nativeImage);
        // vector uVec and vVec define a plane through the origin, 
        // perpendicular to the view vector viewVec
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);
        
        // image is square
        int imageCenter = nativeImage.getWidth() / 2;
        
        // coordinate of point (i,j) in the coordinate system of the volume
        double[] coordOnPlane = new double[3];
        // coordinate of point on the ray casting line
        double[] voxelCoord = new double[3];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);
        
        // sample on a plane through the origin of the volume data
        TFColor voxelColor = new TFColor();

        
        for (int j = 0; j < nativeImage.getHeight(); j++) {
            for (int i = 0; i < nativeImage.getWidth(); i++) {
                int val = 0;
                
                coordOnPlane[0] = uVec[0] * (i - imageCenter) * (1/renderScale) + vVec[0] * (j - imageCenter) * (1/renderScale) + volumeCenter[0];
                coordOnPlane[1] = uVec[1] * (i - imageCenter) * (1/renderScale) + vVec[1] * (j - imageCenter) * (1/renderScale) + volumeCenter[1];
                coordOnPlane[2] = uVec[2] * (i - imageCenter) * (1/renderScale) + vVec[2] * (j - imageCenter) * (1/renderScale) + volumeCenter[2];
                
                double[] fbDepth = getRayDepth(viewVec, coordOnPlane);
                double gap = (fbDepth[0] - fbDepth[1]) / sampleNum;

                // System.out.println("mip::gap::fbDepth " + Arrays.toString(fbDepth));
                // ray cast through (i,j)
                for (double k = fbDepth[1]; k < fbDepth[0]; k += gap) {
                    voxelCoord[0] = coordOnPlane[0] + viewVec[0] * k;
                    voxelCoord[1] = coordOnPlane[1] + viewVec[1] * k;
                    voxelCoord[2] = coordOnPlane[2] + viewVec[2] * k;
                    
                    val = Math.max(val, getVoxel(voxelCoord));
                }
                
                // Map the intensity to a grey value by linear scaling
                //voxelColor.r = val/max;
                //voxelColor.g = voxelColor.r;
                //voxelColor.b = voxelColor.r;
                //voxelColor.a = val > 0 ? 1.0 : 0.0;  // this makes intensity 0 completely transparent and the rest opaque
                // Alternatively, apply the transfer function to obtain a color
                voxelColor = tFunc.getColor(val);
                
                
                // BufferedImage expects a pixel color packed as ARGB in an int
                int c_alpha = voxelColor.a <= 1.0 ? (int) Math.floor(voxelColor.a * 255) : 255;
                int c_red = voxelColor.r <= 1.0 ? (int) Math.floor(voxelColor.r * 255) : 255;
                int c_green = voxelColor.g <= 1.0 ? (int) Math.floor(voxelColor.g * 255) : 255;
                int c_blue = voxelColor.b <= 1.0 ? (int) Math.floor(voxelColor.b * 255) : 255;
                int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
                nativeImage.setRGB(i, j, pixelColor);
            }
        }
        
        scaleImageTo(nativeImage, image);
    }
    
    
    TFColor compositingColor(TFColor oldc, TFColor newc){
        TFColor cColor = new TFColor();
    	cColor.r = newc.a*newc.r + (1-newc.a)*oldc.r;
    	cColor.g = newc.a*newc.g + (1-newc.a)*oldc.g;
    	cColor.b = newc.a*newc.b + (1-newc.a)*oldc.b;
    	return cColor;
    }
    
    void compositing(double[] viewMatrix){
        clearImage(nativeImage);

        // vector uVec and vVec define a plane through the origin, 
        // perpendicular to the view vector viewVec
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);
        
        // image is square
        int imageCenter = nativeImage.getWidth() / 2;
        
        // coordinate of point (i,j) in the coordinate system of the volume
        double[] coordOnPlane = new double[3];
        // coordinate of point on the ray casting line
        double[] voxelCoord = new double[3];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);
        
        // sample on a plane through the origin of the volume data
        double max = volume.getMaximum();

        
        for (int j = 0; j < nativeImage.getHeight(); j++) {
            for (int i = 0; i < nativeImage.getWidth(); i++) {
                int val = 0;
                
                coordOnPlane[0] = uVec[0] * (i - imageCenter) * (1/renderScale) + vVec[0] * (j - imageCenter) * (1/renderScale) + volumeCenter[0];
                coordOnPlane[1] = uVec[1] * (i - imageCenter) * (1/renderScale) + vVec[1] * (j - imageCenter) * (1/renderScale) + volumeCenter[1];
                coordOnPlane[2] = uVec[2] * (i - imageCenter) * (1/renderScale) + vVec[2] * (j - imageCenter) * (1/renderScale)+ volumeCenter[2];
                
                double[] fbDepth = getRayDepth(viewVec, coordOnPlane);
                double gap = (fbDepth[0] - fbDepth[1]) / sampleNum;
                TFColor voxelColor = tFunc.getColor(0);
                // System.out.println("mip::gap: " + gap);
                // ray cast through (i,j)
                for (double k = fbDepth[1]; k < fbDepth[0]; k += gap) {
                    voxelCoord[0] = coordOnPlane[0] + viewVec[0] * k;
                    voxelCoord[1] = coordOnPlane[1] + viewVec[1] * k;
                    voxelCoord[2] = coordOnPlane[2] + viewVec[2] * k;
                    val = getVoxel(voxelCoord);
                    voxelColor = compositingColor(voxelColor, tFunc.getColor(val));
                }
                
                // Map the intensity to a grey value by linear scaling
                // voxelColor.r = val/max;
                // voxelColor.g = voxelColor.r;
                // voxelColor.b = voxelColor.r;
                // voxelColor.a = val > 0 ? 1.0 : 0.0;  // this makes intensity 0 completely transparent and the rest opaque
                // Alternatively, apply the transfer function to obtain a color
                // voxelColor = tFunc.getColor(val);
                
                // BufferedImage expects a pixel color packed as ARGB in an int
                int c_alpha = voxelColor.a <= 1.0 ? (int) Math.floor(voxelColor.a * 255) : 255;
                int c_red = voxelColor.r <= 1.0 ? (int) Math.floor(voxelColor.r * 255) : 255;
                int c_green = voxelColor.g <= 1.0 ? (int) Math.floor(voxelColor.g * 255) : 255;
                int c_blue = voxelColor.b <= 1.0 ? (int) Math.floor(voxelColor.b * 255) : 255;
                int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
                nativeImage.setRGB(i, j, pixelColor);
            }
        }
        
        scaleImageTo(nativeImage, image);
    }
    
    private double levoy(double[] p, int v) {
    	double f = v;
    	double fv = tfEditor2D.triangleWidget.baseIntensity;
    	double g = gradients.getInterpolatedGrad(p[0], p[1], p[2]).mag;
    	if(g > tfEditor2D.triangleWidget.maxGrad || g < tfEditor2D.triangleWidget.minGrad)
    		return 0.0;
    	double r = tfEditor2D.triangleWidget.radius;
    	double alpha;
    	if(g == 0 && f == fv)
    		alpha = 1;
    	else if(g > 0 && f-r*g <= fv && fv <= f+r*g)
    		alpha = 1 - ((1/r)*Math.abs((fv-f)/(g)));
    	else
    		alpha = 0;
    	return (tfEditor2D.triangleWidget.color.a)*alpha;
    }
    
    private TFColor deepCopyTFColor(TFColor source){
        return new TFColor(source.r, source.g, source.b, source.a);
    }
    private TFColor phongShading(double[] L, double[] p){
        // head light assumption H = L = V
        double[] N = new double[3];
        // reverse view vector
        VectorMath.setVector(L, -L[0], -L[1], -L[2]);
        VoxelGradient grad = gradients.getInterpolatedGrad(p[0], p[1], p[2]);
        TFColor c = deepCopyTFColor(tfEditor2D.triangleWidget.color);
        if (grad.mag == 0){
            N[0] = 0;
            N[1] = 0;
            N[2] = 0;
        } else {
            N[0] = grad.x/grad.mag;
            N[1] = grad.y/grad.mag;
            N[2] = grad.z/grad.mag;
        }
        double LN = Math.abs(VectorMath.dotproduct(L, N));
        double NH = Math.abs(VectorMath.dotproduct(N, L));
        
        double ambient = kambient;
        double specular = kspec * Math.pow(NH, phongAlpha);

        c.r = c.r*kdiff*LN + ambient + specular;
        c.g = c.g*kdiff*LN + ambient + specular;
        c.b = c.b*kdiff*LN + ambient + specular;
        // c.a remains
        return c;
    }
    
    void transferFunc2D(double[] viewMatrix){
        clearImage(nativeImage);

        // vector uVec and vVec define a plane through the origin, 
        // perpendicular to the view vector viewVec
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);
        
        // image is square
        int imageCenter = nativeImage.getWidth() / 2;
        
        // coordinate of point (i,j) in the coordinate system of the volume
        double[] coordOnPlane = new double[3];
        // coordinate of point on the ray casting line
        double[] voxelCoord = new double[3];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);
        
        // sample on a plane through the centre of the volume data
        // double max = volume.getMaximum();
        
        for (int j = 0; j < nativeImage.getHeight(); j++) {
            for (int i = 0; i < nativeImage.getWidth(); i++) {
                int val = 0;
                
                coordOnPlane[0] = uVec[0] * (i - imageCenter) * (1/renderScale) + vVec[0] * (j - imageCenter) * (1/renderScale) + volumeCenter[0];
                coordOnPlane[1] = uVec[1] * (i - imageCenter) * (1/renderScale) + vVec[1] * (j - imageCenter) * (1/renderScale) + volumeCenter[1];
                coordOnPlane[2] = uVec[2] * (i - imageCenter) * (1/renderScale) + vVec[2] * (j - imageCenter) * (1/renderScale)+ volumeCenter[2];
                
                double[] fbDepth = getRayDepth(viewVec, coordOnPlane);
                double gap = (fbDepth[0] - fbDepth[1]) / sampleNum;
                TFColor col = deepCopyTFColor(tfEditor2D.triangleWidget.color);
                double alpha = 1;
                // ray cast through (i,j)
                for (double k = fbDepth[1]; k < fbDepth[0]; k += gap) {
                    voxelCoord[0] = coordOnPlane[0] + viewVec[0] * k;
                    voxelCoord[1] = coordOnPlane[1] + viewVec[1] * k;
                    voxelCoord[2] = coordOnPlane[2] + viewVec[2] * k;
                    val = getVoxel(voxelCoord);
                    //grad = gradients.getInterpolatedGrad(voxelCoord[0],voxelCoord[1],voxelCoord[2]).mag;
                    double alpha1 = levoy(voxelCoord, val);
                    alpha = alpha*(1-alpha1);
                }
                alpha = 1 - alpha;
                if (enableShading){
                    col = phongShading(viewVec, coordOnPlane);
                    //col = anotherShading(viewVec, voxelCoord, deepCopyTFColor(col));
                }
                // Map the intensity to a grey value by linear scaling
                // voxelColor.r = val/max;
                // voxelColor.g = voxelColor.r;
                // voxelColor.b = voxelColor.r;
                // voxelColor.a = val > 0 ? 1.0 : 0.0;  // this makes intensity 0 completely transparent and the rest opaque
                // Alternatively, apply the transfer function to obtain a color
                // voxelColor = tFunc.getColor(val);
                
                // BufferedImage expects a pixel color packed as ARGB in an int
                int c_alpha = alpha <= 1.0 ? (int) Math.floor(alpha* 255) : 255;
                if (col.a < 0){
                    System.out.println("alpha wrong");
                    c_alpha = 0;
                }
                int c_red = col.r <= 1.0 ? (int) Math.floor( col.r * 255) : 255;
                int c_green = col.g <= 1.0 ? (int) Math.floor( col.g * 255) : 255;
                int c_blue = col.b <= 1.0 ? (int) Math.floor( col.b * 255) : 255;
                int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
                nativeImage.setRGB(i, j, pixelColor);
            }
        }
        
        scaleImageTo(nativeImage, image);
    }
    
    public void setRayFunction(int functionName){
        this.rayFunction = functionName;
    }
    
    private void callRayFunction() {
        switch(this.rayFunction) {
            case 0:
                slicer(viewMatrix);
                break;
            case 1:
                mip(viewMatrix);
                break;
            case 2:
                compositing(viewMatrix);
                break;
            case 3:
                transferFunc2D(viewMatrix);
                break;
            default:
                slicer(viewMatrix);
                break;
        }
    }
    
    private void drawBoundingBox(GL2 gl) {
        gl.glPushAttrib(GL2.GL_CURRENT_BIT);
        gl.glDisable(GL2.GL_LIGHTING);
        gl.glColor4d(1.0, 1.0, 1.0, 1.0);
        gl.glLineWidth(1.5f);
        gl.glEnable(GL.GL_LINE_SMOOTH);
        gl.glHint(GL.GL_LINE_SMOOTH_HINT, GL.GL_NICEST);
        gl.glEnable(GL.GL_BLEND);
        gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glDisable(GL.GL_LINE_SMOOTH);
        gl.glDisable(GL.GL_BLEND);
        gl.glEnable(GL2.GL_LIGHTING);
        gl.glPopAttrib();

    }
    
    public void setRenderScale(double rs) {
        renderScale = rs;
        createNativeImage();
    }
    
    private void createNativeImage() {
        nativeImage = new BufferedImage(
            (int) Math.floor(image.getWidth() * renderScale),
            (int) Math.floor(image.getHeight() * renderScale),
            BufferedImage.TYPE_INT_ARGB);
    }
    
    public void setShading(boolean shading){
        this.enableShading = shading;
        callRayFunction();
        createNativeImage();
    }
    
    @Override
    public void visualize(GL2 gl) {


        if (volume == null) {
            return;
        }

        drawBoundingBox(gl);

        gl.glGetDoublev(GL2.GL_MODELVIEW_MATRIX, viewMatrix, 0);

        long startTime = System.currentTimeMillis();
        callRayFunction();    
        
        long endTime = System.currentTimeMillis();
        double runningTime = (endTime - startTime);
        panel.setSpeedLabel(Double.toString(runningTime));

        Texture texture = AWTTextureIO.newTexture(gl.getGLProfile(), image, false);

        gl.glPushAttrib(GL2.GL_LIGHTING_BIT);
        gl.glDisable(GL2.GL_LIGHTING);
        gl.glEnable(GL.GL_BLEND);
        gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);

        // draw rendered image as a billboard texture
        texture.enable(gl);
        texture.bind(gl);
        double halfWidth = image.getWidth() / 2.0;
        gl.glPushMatrix();
        gl.glLoadIdentity();
        gl.glBegin(GL2.GL_QUADS);
        gl.glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
        gl.glTexCoord2d(0.0, 0.0);
        gl.glVertex3d(-halfWidth, -halfWidth, 0.0);
        gl.glTexCoord2d(0.0, 1.0);
        gl.glVertex3d(-halfWidth, halfWidth, 0.0);
        gl.glTexCoord2d(1.0, 1.0);
        gl.glVertex3d(halfWidth, halfWidth, 0.0);
        gl.glTexCoord2d(1.0, 0.0);
        gl.glVertex3d(halfWidth, -halfWidth, 0.0);
        gl.glEnd();
        texture.disable(gl);
        texture.destroy(gl);
        gl.glPopMatrix();

        gl.glPopAttrib();


        if (gl.glGetError() > 0) {
            System.out.println("some OpenGL error: " + gl.glGetError());
        }

    }
    private BufferedImage image;
    private BufferedImage nativeImage;
    private double renderScale = 0.5;
    private int sampleNum = 80;
    private double sampleDepth = 0; // Used for the slicer method
    
    private boolean enableShading = false;
    private double kambient = 0.1;
    private double kdiff = 0.7;
    private double kspec = 0.2;
    private int phongAlpha = 10;

    public void setSampleNum(int sn) {
        sampleNum = sn;
        callRayFunction();
        createNativeImage();
    }
    
    public void setSampleDepth(double sd) {
        sampleDepth = sd;
    }
    
    private void clearImage(BufferedImage image) {
        for (int j = 0; j < image.getHeight(); j++) {
             for (int i = 0; i < image.getWidth(); i++) {
                 image.setRGB(i, j, 0);
             }
         }
    }
    
    private void scaleImageTo(BufferedImage img1, BufferedImage img2) {
        AffineTransform at = new AffineTransform();
        at.scale(1/renderScale, 1/renderScale);
        AffineTransformOp scaleOp = new AffineTransformOp(at, AffineTransformOp.TYPE_NEAREST_NEIGHBOR);
        
        scaleOp.filter(img1, img2);
    }
    
    private double[] viewMatrix = new double[4 * 4];

    @Override
    public void changed() {
        callRayFunction();
        for (int i=0; i < listeners.size(); i++) {
            listeners.get(i).changed();
        }
    }
}
