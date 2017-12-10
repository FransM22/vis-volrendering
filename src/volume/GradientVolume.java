/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package volume;

/**
 *
 * @author michel
 */
public class GradientVolume {

    public GradientVolume(Volume vol) {
        volume = vol;
        dimX = vol.getDimX();
        dimY = vol.getDimY();
        dimZ = vol.getDimZ();
        data = new VoxelGradient[dimX * dimY * dimZ];
        maxmag = -1.0;
        compute();
    }

    public VoxelGradient getGradient(int x, int y, int z) {
        if (x < 0 || x >= volume.getDimX() ||
                y < 0 || y >= volume.getDimY() ||
                z < 0 || z >= volume.getDimZ()) {
            return new VoxelGradient(0, 0, 0);
        }
        return data[x + dimX * (y + dimY * z)];
    }
    
    public VoxelGradient getInterpolatedGrad(double x, double y, double z){
        int x0 = (int)Math.floor(x);
        int x1 = x0 + 1;
        int y0 = (int)Math.floor(y);
        int y1 = y0 +1;
        int z0 = (int)Math.floor(z);
        int z1 = z0 + 1;
        
        VoxelGradient g000 = getGradient(x0, y0, z0);
        VoxelGradient g001 = getGradient(x0, y0, z1);
        VoxelGradient g010 = getGradient(x0, y1, z0);
        VoxelGradient g011 = getGradient(x0, y1, z1);
        VoxelGradient g100 = getGradient(x1, y0, z0);
        VoxelGradient g101 = getGradient(x1, y0, z1);
        VoxelGradient g110 = getGradient(x1, y1, z0);
        VoxelGradient g111 = getGradient(x1, y1, z1);
        
        double ratiox = x - x0;
        double ratioy = y - y0;
        double ratioz = z - z0;
        
        VoxelGradient g00 = linearInterpolate(g000, g001, ratioz);
        VoxelGradient g01 = linearInterpolate(g010, g011, ratioz);
        VoxelGradient g10 = linearInterpolate(g100, g101, ratioz);
        VoxelGradient g11 = linearInterpolate(g110, g111, ratioz);
        
        VoxelGradient g0 = linearInterpolate(g00, g01, ratioy);
        VoxelGradient g1 = linearInterpolate(g10, g11, ratioy);
        
        VoxelGradient g = linearInterpolate(g0, g1, ratiox);
        
        return g;
    }
    
    private VoxelGradient linearInterpolate(VoxelGradient g1, VoxelGradient g2, double r){
        float x = (float)(g1.x*(1-r) + g2.x*r);
        float y = (float)(g1.y*(1-r) + g2.y*r);
        float z = (float)(g1.z*(1-r) + g2.z*r);
        
        return new VoxelGradient(x,y,z);
    }
    
    public void setGradient(int x, int y, int z, VoxelGradient value) {
        data[x + dimX * (y + dimY * z)] = value;
    }

    public void setVoxel(int i, VoxelGradient value) {
        data[i] = value;
    }

    public VoxelGradient getVoxel(int i) {
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

    private void compute() {
        for (int x = 0; x < dimX; x++){
            for (int y = 0; y < dimY; y++){
                for (int z = 0; z < dimZ; z++){
                    float gx, gy, gz;
                    
                    if (x == 0 || x == dimX - 1){
                        gx = 0;
                    } else {
                        gx = ((float)volume.getVoxel(x+1, y, z) - (float)volume.getVoxel(x-1, y, z)) / 2;
                    }
                    
                    if (y == 0 || y == dimY - 1){
                        gy = 0;
                    } else {
                        gy = ((float)volume.getVoxel(x, y+1, z) - (float)volume.getVoxel(x, y-1, z)) / 2;
                    }
                    
                    if (z ==0 || z == dimZ - 1){
                        gz = 0;
                    } else {
                        gz = ((float)volume.getVoxel(x, y, z+1) - (float)volume.getVoxel(x, y, z-1)) / 2;
                    }
                    
                    VoxelGradient vg = new VoxelGradient(gx, gy, gz);
                    setGradient(x, y, z, vg);
                    if (maxmag < vg.mag){
                        maxmag = vg.mag;
                    }
                }
            }
        }
                
    }
    
    public double getMaxGradientMagnitude() {
        if (maxmag >= 0) {
            return maxmag;
        } else {
            double magnitude = data[0].mag;
            for (int i=0; i<data.length; i++) {
                magnitude = data[i].mag > magnitude ? data[i].mag : magnitude;
            }   
            maxmag = magnitude;
            return magnitude;
        }
    }
    
    private int dimX, dimY, dimZ;
    private VoxelGradient zero = new VoxelGradient();
    VoxelGradient[] data;
    Volume volume;
    double maxmag;
}
