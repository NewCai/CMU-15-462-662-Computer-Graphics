#include "software_renderer.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#include "triangulation.h"

using namespace std;

namespace CMU462 {

// Implements SoftwareRenderer //

void SoftwareRendererImp::draw_svg(SVG& svg) {
  // set top level transformation
 // transformation = canvas_to_screen;
  transform_stack.push(canvas_to_screen);
  memset(&supersample_target[0], 0, supersample_target.size() * sizeof supersample_target[0]);
  // draw all elements
  for (size_t i = 0; i < svg.elements.size(); ++i) {
    draw_element(svg.elements[i]);
  }

  transformation = transform_stack.top();

  // draw canvas outline
  Vector2D a = transform(Vector2D(0, 0));
  a.x--;
  a.y--;
  Vector2D b = transform(Vector2D(svg.width, 0));
  b.x++;
  b.y--;
  Vector2D c = transform(Vector2D(0, svg.height));
  c.x--;
  c.y++;
  Vector2D d = transform(Vector2D(svg.width, svg.height));
  d.x++;
  d.y++;

  rasterize_line(a.x, a.y, b.x, b.y, Color::Black);
  rasterize_line(a.x, a.y, c.x, c.y, Color::Black);
  rasterize_line(d.x, d.y, b.x, b.y, Color::Black);
  rasterize_line(d.x, d.y, c.x, c.y, Color::Black);

  // resolve and send to render target
  resolve();
  transform_stack.pop();
}

void SoftwareRendererImp::set_sample_rate(size_t sample_rate) {
  // Task 4:
  // You may want to modify this for supersampling support
  this->sample_rate = sample_rate;
  this->sample_w = this->target_w * this->sample_rate;
  this->sample_h = this->target_h * this->sample_rate;
  this->supersample_target.resize(4 * this->sample_w * this->sample_h);
}

void SoftwareRendererImp::set_render_target(unsigned char* render_target,
                                            size_t width, size_t height) {
  // Task 4:
  // You may want to modify this for supersampling support
  this->render_target = render_target;
  this->target_w = width;
  this->target_h = height;
  this->sample_w = this->target_w * this->sample_rate;
  this->sample_h = this->target_h * this->sample_rate;
  this->supersample_target.resize(4 * this->sample_w * this->sample_h);
}

void SoftwareRendererImp::draw_element(SVGElement* element) {
  // Task 5 (part 1):
  // Modify this to implement the transformation stack
  transform_stack.push(transform_stack.top() * element->transform);
  transformation = transform_stack.top();
  
  switch (element->type) {
    case POINT:
      draw_point(static_cast<Point&>(*element));
      break;
    case LINE:
      draw_line(static_cast<Line&>(*element));
      break;
    case POLYLINE:
      draw_polyline(static_cast<Polyline&>(*element));
      break;
    case RECT:
      draw_rect(static_cast<Rect&>(*element));
      break;
    case POLYGON:
      draw_polygon(static_cast<Polygon&>(*element));
      break;
    case ELLIPSE:
      draw_ellipse(static_cast<Ellipse&>(*element));
      break;
    case IMAGE:
      draw_image(static_cast<Image&>(*element));
      break;
    case GROUP:
      draw_group(static_cast<Group&>(*element));
      break;
    default:
      break;
  }

  transform_stack.pop();
}

// Primitive Drawing //

void SoftwareRendererImp::draw_point(Point& point) {
  Vector2D p = transform(point.position);
  rasterize_point(p.x, p.y, point.style.fillColor);
}

void SoftwareRendererImp::draw_line(Line& line) {
  Vector2D p0 = transform(line.from);
  Vector2D p1 = transform(line.to);
  rasterize_line(p0.x, p0.y, p1.x, p1.y, line.style.strokeColor);
}

void SoftwareRendererImp::draw_polyline(Polyline& polyline) {
  Color c = polyline.style.strokeColor;
  float w = polyline.style.strokeWidth;

  if (c.a != 0) {
    int nPoints = polyline.points.size();
    for (int i = 0; i < nPoints - 1; i++) {
      Vector2D p0 = transform(polyline.points[(i + 0) % nPoints]);
      Vector2D p1 = transform(polyline.points[(i + 1) % nPoints]);
      rasterize_line(p0.x, p0.y, p1.x, p1.y, c, w);
    }
  }
}

void SoftwareRendererImp::draw_rect(Rect& rect) {
  Color c;

  // draw as two triangles
  float x = rect.position.x;
  float y = rect.position.y;
  float w = rect.dimension.x;
  float h = rect.dimension.y;

  Vector2D p0 = transform(Vector2D(x, y));
  Vector2D p1 = transform(Vector2D(x + w, y));
  Vector2D p2 = transform(Vector2D(x, y + h));
  Vector2D p3 = transform(Vector2D(x + w, y + h));

  // draw fill
  c = rect.style.fillColor;
  if (c.a != 0) {
    rasterize_triangle(p0.x, p0.y, p1.x, p1.y, p2.x, p2.y, c);
    rasterize_triangle(p2.x, p2.y, p1.x, p1.y, p3.x, p3.y, c);
  }

  // draw outline
  c = rect.style.strokeColor;
  if (c.a != 0) {
    rasterize_line(p0.x, p0.y, p1.x, p1.y, c);
    rasterize_line(p1.x, p1.y, p3.x, p3.y, c);
    rasterize_line(p3.x, p3.y, p2.x, p2.y, c);
    rasterize_line(p2.x, p2.y, p0.x, p0.y, c);
  }
}

void SoftwareRendererImp::draw_polygon(Polygon& polygon) {
  Color c;

  // draw fill
  c = polygon.style.fillColor;
  if (c.a != 0) {
    // triangulate
    vector<Vector2D> triangles;
    triangulate(polygon, triangles);

    // draw as triangles
    for (size_t i = 0; i < triangles.size(); i += 3) {
      Vector2D p0 = transform(triangles[i + 0]);
      Vector2D p1 = transform(triangles[i + 1]);
      Vector2D p2 = transform(triangles[i + 2]);
      rasterize_triangle(p0.x, p0.y, p1.x, p1.y, p2.x, p2.y, c);
    }
  }
  
  // draw outline
  c = polygon.style.strokeColor;
  if (c.a != 0) {
    int nPoints = polygon.points.size();
    for (int i = 0; i < nPoints; i++) {
      Vector2D p0 = transform(polygon.points[(i + 0) % nPoints]);
      Vector2D p1 = transform(polygon.points[(i + 1) % nPoints]);
      rasterize_line(p0.x, p0.y, p1.x, p1.y, c);
    }
  }
}

void SoftwareRendererImp::draw_ellipse(Ellipse& ellipse) {
  // Extra credit
}

void SoftwareRendererImp::draw_image(Image& image) {
  Vector2D p0 = transform(image.position);
  Vector2D p1 = transform(image.position + image.dimension);

  rasterize_image(p0.x, p0.y, p1.x, p1.y, image.tex);
}

void SoftwareRendererImp::draw_group(Group& group) {
  for (size_t i = 0; i < group.elements.size(); ++i) {
    draw_element(group.elements[i]);
  }
}

// Rasterization //

// The input arguments in the rasterization functions
// below are all defined in screen space coordinates

void SoftwareRendererImp::rasterize_point(float x, float y, Color c) {
    // fill in the nearest pixel
  int sx = (int)floor(x);
  int sy = (int)floor(y);

  // check bounds
  if (sx < 0 || sx >= this->target_w) return;
  if (sy < 0 || sy >= this->target_h) return;

  sx *= sample_rate;
  sy *= sample_rate;

  for (int i = sx; i < sx + sample_rate; ++i) {
    for (int j = sy; j < sy + sample_rate; ++j) {
      size_t pxy = 4 * (i + j * sample_w);
      supersample_target[pxy] = (uint8_t)(c.r * 255);
      supersample_target[pxy + 1] = (uint8_t)(c.g * 255);
      supersample_target[pxy + 2] = (uint8_t)(c.b * 255);
      supersample_target[pxy + 3] = (uint8_t)(c.a * 255);
    }
  }

  #if 0
  float step = 1 / sample_rate;
  for (int i = 0; i < sample_rate; ++i) {
    for (int j = 0; j < sample_rate; ++j) {
      fill_sample(x + i * step, y + j * step, color);
    }
  }
  #endif
}

void SoftwareRendererImp::rasterize_line(float x0, float y0, float x1, float y1,
                                         Color color, float width) {
  // Task 2:
  // Implement line rasterization

#if 0
  // An possible implementation of drawing with thickness.
  // TODO: A fast implementaton of fllling rect
  Vector2D s(x1 - x0, y1 - y0);
  s /= s.norm();
  s *= 0.6f;
  Vector2D p(-s.y, s.x);
  float xx1 = x0 + p.x;
  float yy1 = y0 + p.y;
  float xx2 = x1 + p.x;
  float yy2 = y1 + p.y;
  float xx3 = x1 - p.x;
  float yy3 = y1 - p.y;
  float xx4 = x0 - p.x;
  float yy4 = y0 - p.y;

  rasterize_line_xiaolinwu(xx1, yy1, xx2, yy2, color, width);
  rasterize_line_xiaolinwu(xx2, yy2, xx3, yy3, color, width);
  rasterize_line_xiaolinwu(xx3, yy3, xx4, yy4, color, width);
  rasterize_line_xiaolinwu(xx4, yy4, xx1, yy1, color, width);
#endif

  rasterize_line_xiaolinwu(x0, y0, x1, y1, color, width);
}

void SoftwareRendererImp::rasterize_line_bresenham(float x0, float y0, float x1, float y1,
                                                   Color color, float width) {
  
  x0 *= sample_rate;
  x1 *= sample_rate;
  y0 *= sample_rate;
  y1 *= sample_rate;

  bool steep = abs(x1 - x0) < abs(y1 - y0);
  if (steep) {  // abs(slope) > 1
    swap(x0, y0);
    swap(x1, y1);
  }

  if (x0 > x1) {
    swap(x0, x1);
    swap(y0, y1);
  }

  int dx = x1 - x0, dy = y1 - y0, y = y0, eps = 0;
  int sign = dy > 0 ? 1 : -1;
  dy = abs(dy);
  for (int x = x0; x <= x1; x++) {
    if (steep) fill_sample(y, x, color);
    else fill_sample(x, y, color);
    eps += dy;
    if ((eps << 1) >= dx) {
      y += sign;
      eps -= dx;
    }
  }
}

inline float ipart(float x) {
  return floor(x);
}

inline float round(float x) {
  return ipart(x + 0.5f);
}

inline float fpart(float x) {
  return x - floor(x);
}

inline float rfpart(float x) {
  return 1 - fpart(x);
}

void SoftwareRendererImp::rasterize_line_xiaolinwu(float x0, float y0, float x1, float y1,
                                                   Color color, float width) {
  x0 *= sample_rate;
  x1 *= sample_rate;
  y0 *= sample_rate;
  y1 *= sample_rate;
  
  bool steep = abs(x1 - x0) < abs(y1 - y0);
  if (steep) {  // abs(slope) > 1
    swap(x0, y0);
    swap(x1, y1);
  }

  if (x0 > x1) {
    swap(x0, x1);
    swap(y0, y1);
  }

  float dx = x1 - x0;
  float dy = y1 - y0;
  float gradient;
  if (dx == 0.0f)
    gradient = 1.0;
  else
    gradient = dy / dx;

  // handle first endpoint
  float xend = round(x0);
  float yend = y0 + gradient * (xend - x0);
  float xgap = rfpart(x0 + 0.5f);
  float xpxl1 = xend;  // this will be used in the main loop
  float ypxl1 = ipart(yend);
  
  if (steep) {
    color.a = rfpart(yend) * xgap;
    fill_sample(ypxl1, xpxl1, color);
    color.a = fpart(yend) * xgap;
    fill_sample(ypxl1 + 1, xpxl1, color);
  } else {
    color.a = rfpart(yend) * xgap;
    fill_sample(xpxl1, ypxl1, color);
    color.a = fpart(yend) * xgap;
    fill_sample(xpxl1, ypxl1 + 1, color);
  }

  float intery = yend + gradient;  // first y-intersection for the main loop

  // handle first endpoint
  xend = round(x1);
  yend = y1 + gradient * (xend - x1);
  xgap = fpart(x1 + 0.5f);
  float xpxl2 = xend;  // this will be used in the main loop
  float ypxl2 = ipart(yend);
  if (steep) {
    color.a = rfpart(yend) * xgap;
    fill_sample(ypxl2, xpxl2, color);
    color.a = fpart(yend) * xgap;
    fill_sample(ypxl2 + 1, xpxl2, color);
  } else {
    color.a = rfpart(yend) * xgap;
    fill_sample(xpxl2, ypxl2, color);
    color.a = fpart(yend) * xgap;
    fill_sample(xpxl2, ypxl2 + 1, color);
  }
  
  if (steep) {
    for (float x = xpxl1 + 1; x <= xpxl2 - 1 * sample_rate; ++x) {
      color.a = rfpart(intery);
      fill_sample(ipart(intery), x, color);
      color.a = fpart(intery);
      fill_sample(ipart(intery) + 1, x, color);
      intery += gradient;
    }
  } else {
    for (float x = xpxl1 + 1; x <= xpxl2 - 1 * sample_rate; ++x){
      color.a = rfpart(intery);
      fill_sample(x, ipart(intery), color);
      color.a = fpart(intery);
      fill_sample(x, ipart(intery) + 1, color);
      intery += gradient;
    }
  }
}

inline float test_line(float x, float y, float dy, float dx, float yi, float xi) {
  return (x - xi) * dy - (y - yi) * dx; // N * P < 0
}

void SoftwareRendererImp::rasterize_triangle(float x0, float y0, float x1,
                                             float y1, float x2, float y2,
                                             Color color) {
  x0 *= sample_rate;
  y0 *= sample_rate;
  x1 *= sample_rate;
  y1 *= sample_rate;
  x2 *= sample_rate;
  y2 *= sample_rate;
  // Task 3:
  // Implement triangle rasterization
  float dy10 = y1 - y0, dx10 = x1 - x0,
        dy21 = y2 - y1, dx21 = x2 - x1,
        dy02 = y0 - y2, dx02 = x0 - x2;
  
  float endx = ceil(max(max(x0, x1), x2));
  float endy = ceil(max(max(y0, y1), y2));
  float sx = floor(min(min(x0, x1), x2)) + 0.5f;
  float sy = floor(min(min(y0, y1), y2)) + 0.5f;

  for (float y = sy; y <= endy; y += 1) {
    bool isInTriangle = false;
    for (float x = sx; x <= endx; x += 1) {
      float t1 = test_line(x, y, dy10, dx10, y0, x0);
      float t2 = test_line(x, y, dy21, dx21, y1, x1);
      float t3 = test_line(x, y, dy02, dx02, y2, x2);

      float t12 = t1 * t2, t23 = t2 * t3;
      if (t12 < 0 || t23 < 0) {
        if (isInTriangle) {
          break;
        }
      } else {
        isInTriangle = true;
        fill_sample(x, y, color);
      }
    }
  }
}

void SoftwareRendererImp::rasterize_image(float x0, float y0, float x1,
                                          float y1, Texture& tex) {
  // Task 6:
  // Implement image rasterization
  x0 *= sample_rate;
  x1 *= sample_rate;
  y0 *= sample_rate;
  y1 *= sample_rate;

  float w = x1 - x0, h = y1 - y0;
  float su = 1.0f / w, sv = 1.0f / h;
  float u_scale = tex.width / w, v_scale = tex.height / h;

  float d = max(0.f, log2f(max(u_scale, v_scale)));
  int k0 = floor(d), k1 = k0 + 1;

  for (float x = x0, u = 0; x <= x1; ++x, u += su) {
    for (float y = y0, v = 0; y <= y1; ++y, v += sv) {
      fill_sample(x, y, sampler->sample_trilinear(tex, u, v, u_scale, v_scale));
    }
  }
}

// resolve samples to render target
void SoftwareRendererImp::resolve(void) {
  // Task 4:
  // Implement supersampling
  // You may also need to modify other functions marked with "Task 4".
  size_t sampleNum = sample_rate * sample_rate;
  for (size_t x = 0; x <= sample_w - sample_rate; x += sample_rate) {
    for (size_t y = 0; y <= sample_h - sample_rate; y += sample_rate) {
      uint16_t r = 0, g = 0, b = 0, a = 0;
      for (size_t i = 0; i < sample_rate; ++i) {
        for (size_t j = 0; j < sample_rate; ++j) {
          size_t samplePos = 4 * (x + i + (y +  j) * sample_w);
          r += supersample_target[samplePos];
          g += supersample_target[samplePos + 1];
          b += supersample_target[samplePos + 2];
          a += supersample_target[samplePos + 3];
        }
      }

      r /= sampleNum; g /= sampleNum; b /= sampleNum; a /= sampleNum; 
      size_t sx = x / sample_rate;
      size_t sy = y / sample_rate;
      size_t pixPos = 4 * (sx + sy * target_w);
      render_target[pixPos] = (uint8_t)(r);
      render_target[pixPos + 1] = (uint8_t)(g);
      render_target[pixPos + 2] = (uint8_t)(b);
      render_target[pixPos + 3] = (uint8_t)(a);
    }
  }
}

void SoftwareRendererImp::fill_sample( int x, int y, const Color& c ) {
    // check bounds
  if (x < 0 || x >= this->sample_w) return;
  if (x < 0 || y >= this->sample_h) return;

  // fill sample
  size_t pos = 4 * (x + y * sample_w);

  Color from = c;

  Color scr;
  scr.r = supersample_target[pos] / 255.0f;
  scr.g = supersample_target[pos + 1] / 255.0f;
  scr.b = supersample_target[pos + 2] / 255.0f;
  scr.a = supersample_target[pos + 3] / 255.0f;

  Color to;
  to.r = (1.0f - from.a) * scr .r + from.r * from.a;
  to.g = (1.0f - from.a) * scr .g + from.g * from.a;
  to.b = (1.0f - from.a) * scr .b + from.b * from.a;
  to.a = 1.0f - (1.0f  - from.a) * (1.0f - scr.a);

  supersample_target[pos] = (uint8_t)(to.r * 255);
  supersample_target[pos + 1] = (uint8_t)(to.g * 255);
  supersample_target[pos + 2] = (uint8_t)(to.b * 255);
  supersample_target[pos + 3] = (uint8_t)(to.a * 255);
}

}  // namespace CMU462
