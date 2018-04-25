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

float interpolate_linear(const VRVOL* vol, glm::vec3 pt) {
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
  return linear;
}

float triCubic (glm::vec3 p, const VRVOL* vol, int xDim, int yDim, int zDim)
{
  int             x, y, z;
  register int    i, j, k;
  float           dx, dy, dz;
  register int    pv;
  float           u[4], v[4], w[4];
  float           r[4], q[4];
  float           vox = 0;
  int             xyDim;

  xyDim = xDim * yDim;

  // x = (int) p.x, y = (int) p.y, z = (int) p.z;
  x = (int) floor(p.x), y = (int) floor(p.y), z = (int) floor(p.z);
  if (x < 0 || x >= xDim || y < 0 || y >= yDim || z < 0 || z >= zDim)
    return (0);

  dx = p.x - (float) x, dy = p.y - (float) y, dz = p.z - (float) z;

  int index = ((x - 1) + (y - 1) * xDim + (z - 1) * xyDim);
  pv = index;

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
        r[j] += u[i] * ((uint8_t*)vol->data)[pv]/MAX_UI8;
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

float interpolate_cubic(const VRVOL* vol, glm::vec3 pt) {
  float cubic = triCubic(pt, vol, vol->gridx, vol->gridy, vol->gridz);
  return cubic;
}
