#include "triangle.h"

#include "CMU462/CMU462.h"
#include "GL/glew.h"
#include "math.h"

namespace CMU462 {
namespace StaticScene {

Triangle::Triangle(const Mesh* mesh, vector<size_t>& v) : mesh(mesh), v(v) {}
Triangle::Triangle(const Mesh* mesh, size_t v1, size_t v2, size_t v3)
    : mesh(mesh), v1(v1), v2(v2), v3(v3) {}

BBox Triangle::get_bbox() const
{
	// compute the bounding box of the triangle
	BBox box;

	Vector3D points[3] = {
		mesh->positions[v1],
		mesh->positions[v2],
		mesh->positions[v3]};

	for (size_t i = 0; i < 3; i++)
	{
		Vector3D p = points[i];
		for (size_t j = 0; j < 3; j++)
		{
			box.min[j] = min(box.min[j], p[j]);
			box.max[j] = max(box.max[j], p[j]);
		}
	}

	box.extent = box.max - box.min;

	return box;
}

bool Triangle::intersect(const Ray& r) const {
  // implement ray-triangle intersection. 

	auto p0 = mesh->positions[v1];
	auto p1 = mesh->positions[v2];
	auto p2 = mesh->positions[v3];

	auto e1 = p1 - p0;
	auto e2 = p2 - p0;
	
	auto e1_x_d = cross(e1, r.d);
	auto denominator = dot(e1_x_d, e2);

	if (abs(denominator) <= 1e-4)
		return false;

	double denominator_i = 1.0 / denominator;

	auto s = r.o - p0;

	auto s_x_e2 = cross(e2, s);

	auto t = dot(s_x_e2, e1) * denominator_i;
	if (t < r.min_t || t > r.max_t)
		return false;

	auto u = dot(s_x_e2, r.d) * denominator_i;
	if (u < 0 || u > 1)
		return false;

	auto v = dot(e1_x_d, s) * denominator_i;
	if (v < 0 || v > 1 || u + v > 1)
		return false;

	r.max_t = t;

	return true;
}

bool Triangle::intersect(const Ray& r, Intersection* isect) const {
  // implement ray-triangle intersection. When an intersection takes
  // place, the Intersection data should be updated accordingly

  // implement ray-triangle intersection. 
  // When an intersection takes place, 
  // the Intersection data should be updated accordingly

	auto p0 = mesh->positions[v1];
	auto p1 = mesh->positions[v2];
	auto p2 = mesh->positions[v3];

	auto e1 = p1 - p0;
	auto e2 = p2 - p0;

	auto e1_x_d = cross(e1, r.d);
	double denominator = dot(e1_x_d, e2);

	if (abs(denominator) <= 1e-4)
		return false;

	auto denominator_i = 1.0 / denominator;

	auto s = r.o - p0;

	auto s_x_e2 = cross(e2, s);

	auto t = dot(s_x_e2, e1) * denominator_i;
	if (t < r.min_t || t > r.max_t)
		return false;

	auto u = dot(s_x_e2, r.d) * denominator_i;
	if (u < 0 || u > 1)
		return false;

	auto v = dot(e1_x_d, s) * denominator_i;
	if (v < 0 || v > 1 || u + v > 1)
		return false;

	r.max_t = t;

	auto n0 = mesh->normals[v1];
	auto n1 = mesh->normals[v2];
	auto n2 = mesh->normals[v3];
	auto n = (1 - u - v) * n0 + u * n1 + v * n2;

	isect->t = t;
	isect->n = dot(n, r.d) < 0 ? n : -n;
	isect->primitive = this;
	isect->bsdf = mesh->get_bsdf();

	return true;
}

void Triangle::draw(const Color& c) const {
  glColor4f(c.r, c.g, c.b, c.a);
  glBegin(GL_TRIANGLES);
  glVertex3d(mesh->positions[v1].x, mesh->positions[v1].y,
             mesh->positions[v1].z);
  glVertex3d(mesh->positions[v2].x, mesh->positions[v2].y,
             mesh->positions[v2].z);
  glVertex3d(mesh->positions[v3].x, mesh->positions[v3].y,
             mesh->positions[v3].z);
  glEnd();
}

void Triangle::drawOutline(const Color& c) const {
  glColor4f(c.r, c.g, c.b, c.a);
  glBegin(GL_LINE_LOOP);
  glVertex3d(mesh->positions[v1].x, mesh->positions[v1].y,
             mesh->positions[v1].z);
  glVertex3d(mesh->positions[v2].x, mesh->positions[v2].y,
             mesh->positions[v2].z);
  glVertex3d(mesh->positions[v3].x, mesh->positions[v3].y,
             mesh->positions[v3].z);
  glEnd();
}

}  // namespace StaticScene
}  // namespace CMU462
