#include "vrvolume.h"

#define MAX_UI8 255.0

float neighup(float x) {
  return ceilf(x);
}

float neighdown(float x) {
  return floorf(x);
}

/*
 *  Returns the value of the grid point p
 */
float valforqxxx (const VRVOL* vol, glm::vec3 p)
{
    if(p.x < 0 || p.y < 0 || p.z < 0 ||
       p.x >= vol->gridx || p.y >= vol->gridy || p.z >= vol->gridz)
    {
        return 0.0; //this treats everything outside the volume as 0
                    //there are other choices, but this will work for
                    //most cases.
    }

    int ind = (int)(p.x);
    ind += (int)(p.y)*vol->gridx;
    ind += (int)(p.z)*vol->gridx*vol->gridy;
    float val =  ((uint8_t*)vol->data)[ind]/MAX_UI8;
    return val;
}

float lerp(float x, float x1, float x2, float f00, float f01) {
  return ((x2 - x) / (x2 - x1)) * f00 + ((x - x1) / (x2 - x1)) * f01;
}

float biLerp(float x, float y, float q11, float q12, float q21, float q22, float x1, float x2, float y1, float y2) {
  float r1 = lerp(x, x1, x2, q11, q21);
  float r2 = lerp(x, x1, x2, q12, q22);

  return lerp(y, y1, y2, r1, r2);
}

float triLerp(float x, float y, float z, float q000, float q001, float q010, float q011, float q100, float q101, float q110, float q111, float x1, float x2, float y1, float y2, float z1, float z2) {
  float x00 = lerp(x, x1, x2, q000, q100);
  float x10 = lerp(x, x1, x2, q010, q110);
  float x01 = lerp(x, x1, x2, q001, q101);
  float x11 = lerp(x, x1, x2, q011, q111);
  float r0 = lerp(y, y1, y2, x00, x01);
  float r1 = lerp(y, y1, y2, x10, x11);

  return lerp(z, z1, z2, r0, r1);
}

float TriCubic (glm::vec3 p, const VRVOL* volume, int xDim, int yDim, int zDim)
{
  int             x, y, z;
  register int    i, j, k;
  float           dx, dy, dz;
  register float *pv;
  float           u[4], v[4], w[4];
  float           r[4], q[4];
  float           vox = 0;
  int             xyDim;

  xyDim = xDim * yDim;

  x = (int) p.x, y = (int) p.y, z = (int) p.z;
  if (x < 0 || x >= xDim || y < 0 || y >= yDim || z < 0 || z >= zDim)
    return (0);

  dx = p.x - (float) x, dy = p.y - (float) y, dz = p.z - (float) z;
  int index = ((x - 1) + (y - 1) * xDim + (z - 1) * xyDim);
  pv = &((float*)volume->data)[index];

# define CUBE(x)   ((x) * (x) * (x))
# define SQR(x)    ((x) * (x))
/*
 #define DOUBLE(x) ((x) + (x))
 #define HALF(x)   ...
 *
 * may also be used to reduce the number of floating point
 * multiplications. The IEEE standard allows for DOUBLE/HALF
 * operations.
 */

  /* factors for Catmull-Rom interpolation */

  u[0] = -0.5 * CUBE (dx) + SQR (dx) - 0.5 * dx;
  u[1] = 1.5 * CUBE (dx) - 2.5 * SQR (dx) + 1;
  u[2] = -1.5 * CUBE (dx) + 2 * SQR (dx) + 0.5 * dx;
  u[3] = 0.5 * CUBE (dx) - 0.5 * SQR (dx);

  v[0] = -0.5 * CUBE (dy) + SQR (dy) - 0.5 * dy;
  v[1] = 1.5 * CUBE (dy) - 2.5 * SQR (dy) + 1;
  v[2] = -1.5 * CUBE (dy) + 2 * SQR (dy) + 0.5 * dy;
  v[3] = 0.5 * CUBE (dy) - 0.5 * SQR (dy);

  w[0] = -0.5 * CUBE (dz) + SQR (dz) - 0.5 * dz;
  w[1] = 1.5 * CUBE (dz) - 2.5 * SQR (dz) + 1;
  w[2] = -1.5 * CUBE (dz) + 2 * SQR (dz) + 0.5 * dz;
  w[3] = 0.5 * CUBE (dz) - 0.5 * SQR (dz);

  for (k = 0; k < 4; k++)
  {
    q[k] = 0;
    for (j = 0; j < 4; j++)
    {
      r[j] = 0;
      for (i = 0; i < 4; i++)
      {
        r[j] += u[i] * *pv;
        pv++;
      }
      q[k] += v[j] * r[j];
      pv += xDim - 4;
    }
    vox += w[k] * q[k];
    pv += xyDim - 4 * xDim;
  }
  return (vox < 0 ? 0.0 : vox);
}

float interpolate(const VRVOL* vol, glm::vec3 pt) {
  float x1 = neighup(pt.x);
  float x2 = neighdown(pt.x);

  float y1 = neighup(pt.y);
  float y2 = neighdown(pt.y);

  float z1 = neighup(pt.z);
  float z2 = neighdown(pt.z);

  float q000 = valforqxxx(vol, glm::vec3(x1, y1, z1));

  float q001 = valforqxxx(vol, glm::vec3(x1, y1, z2));
  float q011 = valforqxxx(vol, glm::vec3(x1, y2, z2));

  float q010 = valforqxxx(vol, glm::vec3(x1, y2, z1));
  float q110 = valforqxxx(vol, glm::vec3(x2, y2, z1));

  float q100 = valforqxxx(vol, glm::vec3(x2, y1, z1));
  float q101 = valforqxxx(vol, glm::vec3(x2, y1, z2));

  float q111 = valforqxxx(vol, glm::vec3(x2, y2, z2));

  float linear = triLerp(pt.x, pt.y, pt.z, q000, q001, q010, q011, q100, q101, q110, q111, x1, x2, y1, y2, z1, z2);
  // float cubic = TriCubic(pt, vol, vol->gridx, vol->gridy, vol->gridz);
  // return cubic;
  return linear;
}

/*
 *  Returns the value of the nearest grid point to "pt"
 */
float
interpolate_nearest_ui8(const VRVOL* vol, glm::vec3 pt)
{
    if(pt.x < 0 || pt.y < 0 || pt.z < 0 ||
       pt.x >= vol->gridx || pt.y >= vol->gridy || pt.z >= vol->gridz)
    {
        return 0.0; //this treats everything outside the volume as 0
                    //there are other choices, but this will work for
                    //most cases.
    }

    int ind = (int)roundf(pt.x);
    ind += (int)roundf(pt.y)*vol->gridx;
    ind += (int)roundf(pt.z)*vol->gridx*vol->gridy;
    float val =  ((uint8_t*)vol->data)[ind]/MAX_UI8;
    return val;
}

/*
 * Returns the gradient at pt, this is not the most efficient
 * way to implement this function, but it will work.
 */
glm::vec3
gradient_nearest_ui8(const VRVOL*vol, glm::vec3 pt)
{
    //we can do this more effieinctly by explicitly looking up the data, but
    //this is a better representation.
    glm::vec3 grad;
    grad.x = interpolate_nearest_ui8(vol,pt+glm::vec3(1,0,0))-
                interpolate_nearest_ui8(vol,pt-glm::vec3(1,0,0));
    grad.y = interpolate_nearest_ui8(vol,pt+glm::vec3(0,1,0))-
                interpolate_nearest_ui8(vol,pt-glm::vec3(0,1,0));
    grad.z = interpolate_nearest_ui8(vol,pt+glm::vec3(0,0,1))-
                interpolate_nearest_ui8(vol,pt-glm::vec3(0,0,1));
    return grad;
}


int
vrv_initvolume(VRVOL *vol,
               void *data,
               uint32_t grid_x,
               uint32_t grid_y,
               uint32_t grid_z,
               glm::vec3 vdim,
               VINTERP_T interpolation,
               VOL_T data_type)
{
    if(!data)
    {
        EPRINT("ERROR: Null Data passed to vrc_setvolume\n");
        return 0;
    }
    if(vdim.x <= 0|| vdim.y <=0||vdim.z<=0)
    {
        EPRINT("ERROR: Bad Volume Dimensions w: %f h: %f d: %f\n",vdim.x,vdim.y,vdim.z);
        return 0;
    }

    //compute bounding box
    vol->bb_p0 = glm::vec3(0.0,0.0,0.0);
    vol->bb_p1 = vdim;
    if(vdim.x > 1|| vdim.y > 1 || vdim.z > 1)
    {
        EPRINT("WARNING: Bad Volume Dimensions w: %f h: %f d: %f\n",vdim.x,vdim.y,vdim.z);
        vol->bb_p1 = glm::normalize(vol->bb_p1);
        EPRINT("WARNING: Normalizing to w: %f h: %f d: %f\n",vol->bb_p1.x,vol->bb_p1.y,vol->bb_p1.z);
    }

    //center at <0>
    vol->bb_p0 = -vol->bb_p1* .5f;
    vol->bb_p1 = vol->bb_p1* .5f;

    vol->inv_size = 1.0f/(vol->bb_p1-vol->bb_p0);

    //save data information
    vol->data = data;
    vol->gridx = grid_x;
    vol->gridy = grid_y;
    vol->gridz = grid_z;

    vol->type = data_type;
    vol->interp = interpolation;

    return 1;

}

VRVOL* readvol_ui8(const char *fn, size_t gridx, size_t gridy, size_t gridz, glm::vec3 voldim)
{
    FILE *fd;
    fd = fopen(fn,"rb");
    if(!fd){
        EPRINT("Error opening file %s\n",fn);
    }

    uint8_t *in = (uint8_t *)malloc(gridx*gridy*gridz);

    if(fread(in,1,gridx*gridy*gridz,fd) < gridx*gridy*gridz){
        free(in);
        EPRINT("Error reading in volume %s\n",fn);
    }
    fclose(fd);

    VRVOL *vol = (VRVOL*)malloc(sizeof(VRVOL));
    if(!vrv_initvolume(vol,in,gridx,gridy,gridz,voldim,VI_NEAREST,VRV_UINT8))
    {
        EPRINT("Error creating volume structure\n");
        free(in);
        free(vol);
        return NULL;
    }

    return vol;
}


//ADD a new interpolation scheme by adding it to the
//VINTERP_T enum, for example add VI_LINEAR which will
//get the value 1<<8+1.  Then you can add that combination
//to the switch statment i.e.
//case VRV_UINT8 | VI_LINEAR
//
//Note if you want to support other volume types, say float
//you will need to implement an interpolation function for each
//interpolation type i.e.
//VRV_FLOAT32 | VI_NEAREST and VRV_FLOAT32 | VRV_LINEAR ... etc
//Most people will not need to do this, and it is not required for
//the class
float vrv_interpolate(const VRVOL* vol, glm::vec3 pt)
{
    switch(vol->type|vol->interp){
        case VRV_UINT8 | VI_NEAREST:
            return interpolate(vol, pt); // interpolate_nearest_ui8(vol,pt);
        default:
            EPRINT("BAD VOLUME TYPE or INTERPOLATION METHOD\n");
            return 0;
    }
}

glm::vec3 vrv_gradient(const VRVOL* vol, glm::vec3 pt){
    switch(vol->type|vol->interp){
        case VRV_UINT8|VI_NEAREST:
            return gradient_nearest_ui8(vol,pt);
        default:
            EPRINT("BAD VOLUME IN TYPE or INTERPOLATION Combonation in Volume");
            return glm::vec3(0);
    }
}


/****************************************************************************************************
 * The code was developed by Garrett Aldrich for [ECS 277] Advanced Visualization at UC Davis.
 * Bugs and problems :'(
 * If you are in my class, please don't email me.... start a thread on canvas :-)
 * If you aren't in my class I apologize for the bad code, I probably will never fix it!
 *
 * It's free as in beer
 * (free free, do whatever you want with it)
 *
 * If you use a big chunk, please keep this code block somewhere in your code (the end is fine)
 * Or atleats a comment my name and where you got the code.
 *
 * If you found this useful please don't email me, (sorry I ignore way to many already),
 * feel free to send me a thanks on instagram or something like that (I might not read it for a
 * while, but eventually I will)
 *
 ****************************************************************************************************/
