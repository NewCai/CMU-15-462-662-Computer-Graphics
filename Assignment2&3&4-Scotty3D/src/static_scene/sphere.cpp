#include "sphere.h"

#include <cmath>

#include "../bsdf.h"
#include "../misc/sphere_drawing.h"

namespace CMU462 {
namespace StaticScene {

bool Sphere::test(const Ray& r, double& t1, double& t2) const {
  // TODO (PathTracer):
  // Implement ray - sphere intersection test.
  // Return true if there are intersections and writing the
  // smaller of the two intersection times in t1 and the larger in t2.
  // geometric solution
  auto l = o - r.o; 
  auto tca = dot(l, r.d); 
  // if (tca < 0) return false;
  auto d2 = l.norm2() - tca * tca; 
  if (d2 > r2) return false; 
  auto thc = sqrt(r2 - d2);

  t1 = tca - thc; 
  t2 = tca + thc;
  
  if (t1 > t2) std::swap(t1, t2);

  return true;
}

bool Sphere::intersect(const Ray& r) const {
  // Implement ray - sphere intersection.
  // Note that you might want to use the the Sphere::test helper here.

  double t1, t2;
  if (!test(r, t1, t2))
  {
    return false;
  }

  if (t1 < r.min_t || t2 > r.max_t)
  {
    return false;
  }

  if (t1 < r.max_t)
  {
    r.max_t = t1;
  }
  else
  {
    r.max_t = t2;
  }

  return true;
}

bool Sphere::intersect(const Ray& r, Intersection* isect) const {
  // Implement ray - sphere intersection.
  // Note again that you might want to use the the Sphere::test helper here.
  // When an intersection takes place, the Intersection data should be updated
  // correspondingly.

  double t1, t2;
  if (!test(r, t1, t2))
  {
    return false;
  }

  if (t1 < r.min_t || t2 > r.max_t)
  {
    return false;
  }

  if (t1 < r.max_t)
  {
    r.max_t = t1;
  }
  else
  {
    r.max_t = t2;
  }

	isect->t = r.max_t;
	isect->primitive = this;
	isect->bsdf = object->get_bsdf();
	isect->n = (r.o + r.max_t * r.d - o) / this->r;

	return true;
}

void Sphere::draw(const Color& c) const { Misc::draw_sphere_opengl(o, r, c); }

void Sphere::drawOutline(const Color& c) const {
  // Misc::draw_sphere_opengl(o, r, c);
}

}  // namespace StaticScene
}  // namespace CMU462
