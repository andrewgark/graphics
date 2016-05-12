#ifndef _DRAW
#define _DRAW

#include "GeomSpace.h"
#include <thread>

const int THREADS = 8;
const int ITERATIONS = 4;


std::pair<Point, Figure*> getClosestPoint(GeomSpace &geomSpace, Ray &R) {
    return geomSpace.kdTree.getClosestPointOnRay(R);
}

Color getColorRay(GeomSpace &geomSpace, Ray &R, int iter) {
    std::pair<Point, Figure*> closest = getClosestPoint(geomSpace, R);
    Point P = closest.first;
    Figure* F = closest.second;
    if (!closest.second)
        return Color();

    long double energy = 0;
    Point norm = F->getNormal(P);
    norm = norm / abs(norm);
    if (norm * (R.b - R.a) > EPS)
        norm = norm * (-1);

    for (int i = 0; i < geomSpace.lights.size(); ++i) {
        Light L = geomSpace.lights[i];
        Ray lightRay(L.p, P);
        if (abs(getClosestPoint(geomSpace, lightRay).first - P) < EPS) {
            energy += std::max((long double)0.0, L.e * cosBetween(norm, lightRay.a - lightRay.b) / dist2(L.p, P));
        }
    }
    Color newColor = F->getColor(P) * energy;
    if (iter < ITERATIONS) {
        Point E = projection(R.a, Line(P, P + norm));
        Ray newR(P, E + (E - R.a));

        return newColor * F->alpha + getColorRay(geomSpace, newR, iter + 1) * (1 - F->alpha);
    }
    return newColor;
}

Color getColorPixel(GeomSpace &geomSpace, int x, int y, int w, int h) {
    Point O = geomSpace.screen->a;
    Point OX = geomSpace.screen->b;
    Point OY = geomSpace.screen->d;
    Point P = O + (OX - O) * ((long double) x / (long double) w) + (OY - O) * ((long double) y / (long double) h);
    Ray R(geomSpace.camera, P);
    return getColorRay(geomSpace, R, 0);
}


void draw(GeomSpace *geomSpace, int minWidth, int minHeight, int maxWidth, int maxHeight, png::image< png::rgb_pixel >* image) {

    for (int x = minWidth; x < maxWidth; ++x) {
        for (int y = minHeight; y < maxHeight; ++y) {
            Color color = getColorPixel(*geomSpace, x, y, geomSpace->width, geomSpace->height);
            (*image)[y][x] = png::rgb_pixel(color.r, color.g, color.b); //png::rgb_pixel(x == 0 ? 255 : 0, y == 0 ? 255 : 0, 0);
        }
    }

}

void drawThePicture(GeomSpace &geomSpace) {
    png::image< png::rgb_pixel > image(geomSpace.width, geomSpace.height);
    std::vector<std::thread> v;

    for (int i = 0; i < THREADS; ++i)
        v.push_back(std::thread(draw, &geomSpace,
                                (i * geomSpace.width) / THREADS, 0,
                                ((i + 1) * geomSpace.width) / THREADS, geomSpace.height,
                                &image));
    for (int i = 0; i < THREADS; ++i)
        v[i].join();
    image.write("rgb.png");
}

#endif
