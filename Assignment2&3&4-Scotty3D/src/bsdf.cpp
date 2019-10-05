#include "bsdf.h"

#include <algorithm>
#include <iostream>
#include <utility>


using std::min;
using std::max;
using std::swap;

namespace CMU462 {

void make_coord_space(Matrix3x3& o2w, const Vector3D& n) {
  Vector3D z = Vector3D(n.x, n.y, n.z);
  Vector3D h = z;
  if (fabs(h.x) <= fabs(h.y) && fabs(h.x) <= fabs(h.z))
    h.x = 1.0;
  else if (fabs(h.y) <= fabs(h.x) && fabs(h.y) <= fabs(h.z))
    h.y = 1.0;
  else
    h.z = 1.0;

  z.normalize();
  Vector3D y = cross(h, z);
  y.normalize();
  Vector3D x = cross(z, y);
  x.normalize();

  o2w[0] = x;
  o2w[1] = y;
  o2w[2] = z;
}

// Diffuse BSDF //

Spectrum DiffuseBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return albedo * (1.0 / PI);
}

Spectrum DiffuseBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
	auto xi1 = (double)rand() / RAND_MAX;
	auto xi2 = (double)rand() / RAND_MAX;

	auto a = sqrt(1 - xi1 * xi1);
	wi->x = a * cos(2 * PI * xi2);
	wi->y = a * sin(2 * PI * xi2);
	wi->z = xi1;
	*pdf = 1 / (2 * PI);

	return albedo * (1.0 / PI);
}

// Mirror BSDF //

Spectrum MirrorBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum MirrorBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
  // Implement MirrorBSDF
  reflect(wo, wi);
  *pdf = 1.0f;
  return 1.0 / abs(wi->z) * reflectance;
}

// Glossy BSDF //

/*
Spectrum GlossyBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum GlossyBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
  *pdf = 1.0f;
  return reflect(wo, wi, reflectance);
}
*/

// Refraction BSDF //

Spectrum RefractionBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum RefractionBSDF::sample_f(const Vector3D& wo, Vector3D* wi,
                                  float* pdf) {
  // TODO (PathTracer):
  // Implement RefractionBSDF
  	*pdf = 1.0f;

	if (!refract(wo, wi, ior))
  {
		reflect(wo, wi);
    return Spectrum();
  }

	return 1.0 / abs(wi->z) * transmittance;
}

// Glass BSDF //

Spectrum GlassBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum GlassBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
  // TODO (PathTracer):
  // Compute fresnel coefficient and either reflect or refract based on it.

	if (!refract(wo, wi, ior)) {
		*pdf = 1.0f;
		reflect(wo, wi);
		return reflectance;
	}

  float eta_i = 1.0f;
  float eta_t = ior;

  if (wo.z < 0)
  {
    swap(eta_i, eta_t);
  }

  float cos_theta_i = std::abs(wo.z);
  float cos_theta_t = std::abs(wi->z);

  float r_parl = ((eta_t * cos_theta_i) - (eta_i * cos_theta_t)) /
                ((eta_t * cos_theta_i) + (eta_i * cos_theta_t));
  float r_perp = ((eta_i * cos_theta_i) - (eta_t * cos_theta_t)) /
                ((eta_i * cos_theta_i) + (eta_t * cos_theta_t));

  float fr = (r_parl * r_parl + r_perp * r_perp) / 2;

  if (((double)rand() / RAND_MAX) < fr) {
		*pdf = fr;
		reflect(wo, wi);
		return 1 / cos_theta_t * fr * reflectance;
	}

	*pdf = 1 - fr;
	return 1 / cos_theta_t * (1 - fr) * transmittance;
}

void BSDF::reflect(const Vector3D& wo, Vector3D* wi) {
  // TODO (PathTracer):
  // Implement reflection of wo about normal (0,0,1) and store result in wi.
  wi->x = -wo.x;
  wi->y = -wo.y;
  wi->z = wo.z;
}

bool BSDF::refract(const Vector3D& wo, Vector3D* wi, float ior) {
  // TODO (PathTracer):
  // Use Snell's Law to refract wo surface and store result ray in wi.
  // Return false if refraction does not occur due to total internal reflection
  // and true otherwise. When dot(wo,n) is positive, then wo corresponds to a
  // ray entering the surface through vacuum.
	float ratio = wo.z > 0 ? 1.0 / ior : ior;

	double discriminant = 1 - (1 - wo.z * wo.z) * ratio * ratio;
	if (discriminant < 0)
		return false;

  auto theta = acos(sqrt(discriminant));

	wi->x = - wo.x;
	wi->y = - wo.y;
	wi->z = (wo.z >= 0 ? -1 : 1) * sqrt(wo.x * wo.x + wo.y * wo.y) / abs(tan(theta));
	wi->normalize();

  return true;
}

// Emission BSDF //

Spectrum EmissionBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum EmissionBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
  *wi = sampler.get_sample(pdf);
  return Spectrum();
}

}  // namespace CMU462
