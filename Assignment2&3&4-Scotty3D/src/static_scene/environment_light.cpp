#include "environment_light.h"
#include <algorithm>    // std::min

namespace CMU462 {
namespace StaticScene {

using namespace std;

EnvironmentLight::EnvironmentLight(const HDRImageBuffer* envMap)
    : envMap(envMap) {
  // TODO: (PathTracer) initialize things here as needed
	auto h = envMap->h;
	auto w = envMap->w;
	pMap = vector<vector<double>>(h, vector<double>(w));
	thetaMap = vector<double>(h);

	auto thetaStep = PI / envMap->h;
	double total = 0;
	double theta = 0;
	for (size_t y = 0; y < h; y++, theta += thetaStep) {
		double rowSum = 0;
		for (size_t x = 0; x < w; x++) {
			auto index = x + y * w;
			pMap[y][x] = envMap->data[index].illum() * sin(theta);
			rowSum += pMap[y][x];
		}
		thetaMap[y] = rowSum;
		total += rowSum;
	}

	for (auto &row : pMap)
	{
		for (auto &p : row)
		{
			p /= total;
		}
	}
	for (auto &row : pMap)
	{
		for (size_t i = 1; i < row.size(); ++i)
		{
			row[i] += row[i - 1];
		}
	}

	for (auto &p : thetaMap)
	{
		p /= total;
	}
	for (size_t i = 1; i < thetaMap.size(); ++i)
	{
		thetaMap[i] += thetaMap[i - 1];
	}
}

Spectrum EnvironmentLight::sample_L(const Vector3D& p, Vector3D* wi,
                                    float* distToLight, float* pdf) const {
  /*
  {
  	double Xi1 = (double)(std::rand()) / RAND_MAX;
  	double Xi2 = (double)(std::rand()) / RAND_MAX;
  
  	double theta = acos(1 - 2 * Xi1);
  	double phi = 2.0 * PI * Xi2;
  
  	wi->x = sinf(theta) * cosf(phi);
  	wi->y = sinf(theta) * sinf(phi);
  	wi->z = cosf(theta);
  }
  */
  
  // importance sampling
  {
	double rowPdf = (double)(std::rand()) / RAND_MAX;
	auto hIndex = std::lower_bound(thetaMap.begin(), thetaMap.end(), rowPdf) - thetaMap.begin();

	auto phiPdf = ((double)(std::rand()) / RAND_MAX) * thetaMap[hIndex];
	auto row = pMap[hIndex];
	auto wIndex = std::lower_bound(row.begin(), row.end(), phiPdf) - row.begin();

	double theta = PI * (double)hIndex / envMap->h;
	double phi = 2.0 * PI * (double)wIndex / envMap->w;

	wi->x = sinf(theta) * cosf(phi);
	wi->y = sinf(theta) * sinf(phi);
	wi->z = cosf(theta);
  }


  *distToLight = FLT_MAX;
  *pdf = 0.25f / PI;

  return sample_dir(Ray(p, *wi));
}

Spectrum EnvironmentLight::sample_dir(const Ray &r) const {
	size_t w = envMap->w;
	size_t h = envMap->h;

	double theta = acos(r.d.y);
	double phi = atan2(r.d.x, -r.d.z) + PI;

	double Xi2 = phi / (2.0 * PI);
	double Xi1 = theta / PI;

	auto u = (envMap->w - 1) * Xi2;
	auto v = (envMap->h - 1) * Xi1;

	size_t x = size_t(u);
	size_t y = size_t(v);

	float u_ratio = u - x;
	float v_ratio = v - y;
	float u_opposite = 1 - u_ratio;
	float v_opposite = 1 - v_ratio;

	auto c1 = envMap->data[x + envMap->w * y];
	auto c2 = envMap->data[x + 1 + envMap->w * y];
	auto c3 = envMap->data[x + envMap->w * (y + 1)];
	auto c4 = envMap->data[x + 1 + envMap->w * (y + 1)];

	return (c1 * u_opposite + c2 * u_ratio) * v_opposite + (c3 * u_opposite + c4 * u_ratio) * v_ratio;
}

}  // namespace StaticScene
}  // namespace CMU462
