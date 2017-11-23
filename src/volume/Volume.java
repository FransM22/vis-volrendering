/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package volume;

import java.io.File;
import java.io.IOException;
import static java.lang.Math.floor;

/**
 *
 * @author michel
 */
public class Volume {
    
    public Volume(int xd, int yd, int zd) {
        data = new short[xd*yd*zd];
        dimX = xd;
        dimY = yd;
        dimZ = zd;
    }
    
    public Volume(File file) {
        
        try {
            VolumeIO reader = new VolumeIO(file);
            dimX = reader.getXDim();
            dimY = reader.getYDim();
            dimZ = reader.getZDim();
            data = reader.getData().clone();
            computeHistogram();
        } catch (IOException ex) {
            System.out.println("IO exception");
        }
        
    }
    
    
    public short getVoxel(int x, int y, int z) {
        return data[x + dimX*(y + dimY * z)];
    }
    
    public void setVoxel(int x, int y, int z, short value) {
        data[x + dimX*(y + dimY*z)] = value;
    }

    public void setVoxel(int i, short value) {
        data[i] = value;
    }
    
    public short getVoxel(int i) {
        return data[i];
    }
    
    public int getDimX() {
        return dimX;
    }
    
    public int getDimY() {
        return dimY;
    }
    
    public int getDimZ() {
        return dimZ;
    }

    public short getMinimum() {
        short minimum = data[0];
        for (int i=0; i<data.length; i++) {
            minimum = data[i] < minimum ? data[i] : minimum;
        }
        return minimum;
    }

    public short getMaximum() {
        short maximum = data[0];
        for (int i=0; i<data.length; i++) {
            maximum = data[i] > maximum ? data[i] : maximum;
        }
        return maximum;
    }
 
    public int[] getHistogram() {
        return histogram;
    }
    
    public double getInterpolated(double x, double y, double z) {
        int x_0 = (int) floor(x);
        int x_1 = x_0 + 1;
        double alpha = x - x_0;
        
        int y_0 = (int) floor(y);
        int y_1 = y_0 + 1;
        double beta = y - y_0;

        int z_0 = (int) floor(z);
        int z_1 = z_0 + 1;
        double gamma = z - z_0;
        
        if (x_1 >= getDimX() || y_1 >= getDimY() || z_1 >= getDimZ()) {
            // Volume has no voxel at this position
            return 0;
        }
        
        double s_x0 = getVoxel(x_0, y_0, z_0);
        double s_x1 = getVoxel(x_1, y_0, z_0);
        double s_x2 = getVoxel(x_0, y_1, z_0);
        double s_x3 = getVoxel(x_1, y_1, z_0);
        
        double s_xyz0 =
                (1 - alpha) * (1 - beta) * s_x0 +
                alpha * (1 - beta) * s_x1 +
                (1 - alpha) * beta * s_x2 +
                alpha*beta*s_x3;
        
        // Note, the s_x* variables are being redefined
        s_x0 = getVoxel(x_0, y_0, z_1);
        s_x1 = getVoxel(x_1, y_0, z_1);
        s_x2 = getVoxel(x_0, y_1, z_1);
        s_x3 = getVoxel(x_1, y_1, z_1);
        
        double s_xyz1 =
                (1 - alpha) * (1 - beta) * s_x0 +
                alpha * (1 - beta) * s_x1 +
                (1 - alpha) * beta * s_x2 +
                alpha*beta*s_x3;
                
        double s_xyz =
                (1-gamma) * s_xyz0 +
                gamma * s_xyz1;
        
        return s_xyz;
    }
    
    private void computeHistogram() {
        histogram = new int[getMaximum() + 1];
        for (int i=0; i<data.length; i++) {
            histogram[data[i]]++;
        }
    }
    
    private int dimX, dimY, dimZ;
    private short[] data;
    private int[] histogram;
}
