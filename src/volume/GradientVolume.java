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
        return data[x + dimX * (y + dimY * z)];
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

        // this just initializes all gradients to the vector (0,0,0)
        //for (int i=0; i<data.length; i++) {
        //    data[i] = zero;
        //}
        
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
