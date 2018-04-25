# Volume Ray Casting

> Volume ray casting, sometimes called volumetric ray casting, volumetric ray tracing, or volume ray marching, is an image-based volume rendering technique. It computes 2D images from 3D volumetric data sets (3D scalar fields). Volume ray casting, which processes volume data, must not be mistaken with ray casting in the sense used in ray tracing, which processes surface data. In the volumetric variant, the computation doesn't stop at the surface but "pushes through" the object, sampling the object along the ray. Unlike ray tracing, volume ray casting does not spawn secondary rays. When the context/application is clear, some authors simply call it ray casting. Because raymarching does not necessarily require an exact solution to ray intersection and collisions, it is suitable for real time computing for many applications for which ray tracing is unsuitable.

> Ray casting. For each pixel of the final image, a ray of sight is shot ("cast") through the volume. At this stage it is useful to consider the volume being touched and enclosed within a bounding primitive, a simple geometric object — usually a cuboid — that is used to intersect the ray of sight and the volume.

- [Wikipedia](https://en.wikipedia.org/wiki/Volume_ray_casting)

## Skeleton and Acknowledgement

The volume renderer skeleton is part of the course [ECS 277] Advanced Visualization at UC Davis, developed by
**Garrett Aldrich**. My contribution is the interpolation schemes and enhancements over it.

## Compile and Execute

In the `src` folder:
```
make
./raycaster ../volumes/fuel_8_64.raw 64 64 64
```

## Nearest, Linear & Cubic

![Nearest](https://github.com/neoblizz/raycaster/blob/master/images/fuel_sample.png)
![Linear](https://github.com/neoblizz/raycaster/blob/master/images/fuel_linear_interp.png)
![Cubic](https://github.com/neoblizz/raycaster/blob/master/images/fuel_cubic.png)
