#ifndef _GEOMSPACE
#define _GEOMSPACE

#include "png++.h"
#include <iostream>
#include <cmath>
#include <fstream>
#include <cassert>
#include <iomanip>
#include <limits>
#include <vector>
#include <algorithm>
#include <math.h>

const long double INF = 1e18;
const long double EPS = 1e-7;
const long double PI = atan(1.0) * 4;

long double abs(long double a) {
    if (a < 0)
        return -a;
    else
        return a;
}

struct Color {
    int r, g, b;

    Color() : r(0), g(0), b(0) {}
    Color(int _r, int _g, int _b) : r(_r), g(_g), b(_b) {}

    long double norm(long double c) {
        return std::max((long double)0, std::min((long double)255, c));
    }

    Color operator* (long double d) {
        return Color(norm(r * d), norm(g * d), norm(b * d));
    }

    Color operator+ (Color c) {
        return Color(norm(r + c.r), norm(g + c.g), norm(b + c.g));
    }

    friend std::ostream &operator<< (std::ostream &output, const Color &P)
    {
        output << P.r << " " << P.g << " " << P.b;
        return output;
    }

    friend std::istream &operator>>(std::istream &input, Color &P)
    {
        input >> P.r >> P.g >> P.b;
        return input;
    }
};

struct Texture {
    png::image< png::rgb_pixel > image;

    Texture() {}
    friend std::istream &operator>>(std::istream &input, Texture &T) {
        std::string fileName;
        input >> fileName;
        T.image = png::image< png::rgb_pixel >(fileName.c_str());
        return input;
    }

    Color getColor(long double x, long double y) {
        //std::cout << x << " " << y << "\n";
        int xInt = ((int) x) % image.get_width();
        int yInt = ((int) y) % image.get_height();
        while(xInt < 0)
            xInt += image.get_width();
        while(yInt < 0)
            yInt += image.get_height();
        //std::cout << xInt << " " << yInt << "\n";
        //std::cout << "end\n";
        return Color(image[xInt][yInt].red, image[xInt][yInt].green, image[xInt][yInt].blue);
    }
};

struct Point {
    long double x, y, z;
    Point() : x(0), y(0), z(0) {}
    Point(long double _x, long double _y, long double _z) : x(_x), y(_y), z(_z) {}

    friend std::ostream &operator<< (std::ostream &output, const Point &P)
    {
        output << P.x << " " << P.y << " " << P.z;
        return output;
    }

    friend std::istream &operator>>(std::istream &input, Point &P)
    {
        input >> P.x >> P.y >> P.z;
        return input;
    }

    Point operator+ (Point B) {
        return Point(x + B.x, y + B.y, z + B.z);
    }

    Point operator- (Point B) {
        return Point(x - B.x, y - B.y, z - B.z);
    }

    Point operator* (long double d) {
        return Point(x * d, y * d, z * d);
    }

    Point operator/ (long double d) {
        return Point(x / d, y / d, z / d);
    }

    long double operator* (Point B) {
        return x * B.x + y * B.y + z * B.z;
    }

    Point operator% (Point B) {
        return Point(y * B.z - z * B.y, z * B.x - x * B.z, x * B.y - y * B.x);
    }
};

long double abs(Point A) {
    return sqrt(A * A);
}

long double dist(Point A, Point B) {
    return abs(B - A);
}

long double dist2(Point A, Point B) {
    return (B - A) * (B - A);
}

struct Ray {
    Point a, b;
    Ray() {}
    Ray(Point _a, Point _b) : a(_a), b(_b) {}
};

struct Line {
    Point a, b;
    Line() {}
    Line(Point _a, Point _b) : a(_a), b(_b) {}
};

struct Intersection {
    bool hasIntersection;
    Point point;
    Intersection() : hasIntersection(false) {}
    Intersection(Point _point) : hasIntersection(true), point(_point) {}
};

long double cosBetween(Point X, Point Y) {
    return X * Y / (abs(X) * abs(Y));
}

bool oneDirection(Point X, Point Y) {
    return (abs(X - Y * (abs(X) / abs(Y))) < EPS);
}

bool onOneSide(Point A, Point B, Point C, Point D) {
    Point X = (C - A) % (B - A);
    Point Y = (D - A) % (B - A);
    return oneDirection(X, Y);
}

bool onLine(Point X, Line L) {
    if (abs(X - L.a) < EPS || abs(X - L.b) < EPS)
        return true;
    return oneDirection(L.b - L.a, X - L.a) || oneDirection(L.a - L.b, X - L.a);
}

bool onSameLine(Point X, Point Y, Point Z) {
    return abs(X - Y) < EPS || onLine(Z, Line(X, Y));
}

Point projection(Point X, Line L) {
    if (onLine(X, L))
        return X;
    Point N = (L.b - L.a) % (X - L.a);
    Point V = (L.b - L.a) % N;
    long double h = abs(N) / abs(L.b - L.a);
    Point W = V * (h / abs(V));

    if (onLine(X + W, L))
        return X + W;
    if (onLine(X - W, L))
        return X - W;

    std::cout << X << "\n";
    std::cout << X + W << "\n";
    std::cout << X - W << "\n";
    std::cout << L.a << "\n";
    std::cout << L.b << "\n\n";
    assert(false);

}

struct Figure {
    Color color;
    long double alpha;
    Texture* texture;
    ~Figure() {}
    Figure() : texture(nullptr) {}
    Figure(long double _alpha, Color _clr) : alpha(_alpha), color(_clr) {}
    Figure(long double _alpha, Texture* _txt) : alpha(_alpha), texture(_txt) {}
    virtual Intersection intersectWithRay(Ray &R) = 0;
    virtual Point getNormal(Point &X) = 0;
    virtual Color getColor(Point &X) = 0;
    virtual Point getCenter() = 0;
    virtual Point getMinBBox() = 0;
    virtual Point getMaxBBox() = 0;
};

struct Triangle : public Figure {
    Point a, b, c;
    Triangle(Point _a, Point _b, Point _c) : a(_a), b(_b), c(_c) {}
    Triangle(Point _a, Point _b, Point _c, long double _alpha, Color _clr) : Figure(_alpha, _clr), a(_a), b(_b), c(_c) {}
    Triangle(Point _a, Point _b, Point _c, long double _alpha, Texture* _txt) : Figure(_alpha, _txt), a(_a), b(_b), c(_c) {}

    Intersection intersectWithRay(Ray &R);
    Point getNormal(Point &X);
    Color getColor(Point &X);
    Point getCenter();
    Point getMinBBox() {
        return Point(std::min(a.x, std::min(b.x, c.x)),
                     std::min(a.y, std::min(b.y, c.y)),
                     std::min(a.z, std::min(b.z, c.z)));
    }
    Point getMaxBBox() {
        return Point(std::max(a.x, std::max(b.x, c.x)),
                     std::max(a.y, std::max(b.y, c.y)),
                     std::max(a.z, std::max(b.z, c.z)));
    }
};

struct Quadrangle : public Figure {
    Point a, b, c, d;
    //Quadrangle() {}
    Quadrangle(Point _a, Point _b, Point _c, Point _d) : a(_a), b(_b), c(_c), d(_d) {}
    Quadrangle(Point _a, Point _b, Point _c, Point _d, long double _alpha, Color _clr) :
        Figure(_alpha, _clr), a(_a), b(_b), c(_c), d(_d) {}
    Quadrangle(Point _a, Point _b, Point _c, Point _d, long double _alpha, Texture* _txt) :
        Figure(_alpha, _txt), a(_a), b(_b), c(_c), d(_d) {}
    Quadrangle(Point P, long double w, long double h, Point camera) {
        Point E((long double)0.0, (long double)0.0, (long double)1.0);
        //std::cout << "E = " << E << "\n";
        Point V = E % (P - camera);
        //std::cout << "V = " << V << "\n";
        V = V * ((w / 2) / abs(V));
        //std::cout << "V = " << V << "\n";
        Point G = V % (P - camera);
        G = G * ((h / 2) / abs(G));
        //std::cout << "G = " << G << "\n";

        a = P + V - G;
        b = P - V - G;
        c = P - V + G;
        d = P + V + G;
    }

    Intersection intersectWithRay(Ray &R);
    Point getNormal(Point &X);
    Color getColor(Point &X);
    Point getCenter();
    Point getMinBBox() {
        return Point(std::min(std::min(a.x, b.x), std::min(c.x, d.x)),
                     std::min(std::min(a.y, b.y), std::min(c.y, d.y)),
                     std::min(std::min(a.z, b.z), std::min(c.z, d.z)));
    }
    Point getMaxBBox() {
        return Point(std::max(std::max(a.x, b.x), std::max(c.x, d.x)),
                     std::max(std::max(a.y, b.y), std::max(c.y, d.y)),
                     std::max(std::max(a.z, b.z), std::max(c.z, d.z)));
    }
};

struct Sphere : public Figure {
    Point c;
    long double r;
    Sphere(Point _c, long double _r) : c(_c), r(_r) {}
    Sphere(Point _c, long double _r, long double _alpha, Color _clr) : Figure(_alpha, _clr), c(_c), r(_r) {}

    Intersection intersectWithRay(Ray &R);
    Point getNormal(Point &X);
    Color getColor(Point &X);
    Point getCenter();
    Point getMinBBox() {
        return Point(c.x - r,
                     c.y - r,
                     c.z - r);
    }
    Point getMaxBBox() {
        return Point(c.x + r,
                     c.y + r,
                     c.z + r);
    }
};

struct Plane {
    Point a, b, c;
    Point n, x;
    long double d;
    Plane(Point _a, Point _b, Point _c) : a(_a), b(_b), c(_c) {
        n = (_b - _a) % (_c - _a);
        x = a;
        d = -(n * x);
    }

    Plane(Triangle F) : Plane(F.a, F.b, F.c) {}
    Plane(Triangle *F) : Plane(F->a, F->b, F->c) {}
};

Point Triangle::getNormal(Point &X) {
    Plane G(this);
    return G.n;
}

Color Triangle::getColor(Point &X) {
    if (!texture)
        return color;
    Point u = (b - a);
    u = u / abs(u);
    Point v = getNormal(a) % u;
    v = v / abs(v);
    return texture->getColor(X * u, X * v);
}

Color Quadrangle::getColor(Point &X) {
    if (!texture)
        return color;
    Point u = (b - a);
    u = u / abs(u);
    Point v = getNormal(a) % u;
    v = v / abs(v);
    return texture->getColor(X * u, X * v);
}

Color Sphere::getColor(Point &X) {
    return color;
}


Point Triangle::getCenter() {
    return (a + b + c) / 3;
}

Point Quadrangle::getNormal(Point &X) {
    Triangle *T = new Triangle(a, b, c);
    Plane G(T);
    delete(T);
    return G.n;
}

Point Quadrangle::getCenter() {
    return (a + b + c + d) / 4;
}

Point Sphere::getNormal(Point &X) {
    return X - c;
}

Point Sphere::getCenter() {
    return c;
}

bool isInsideTriangle(Point &P, Triangle F) {
    return (onOneSide(F.a, F.b, P, F.c) && onOneSide(F.a, F.c, P, F.b) && onOneSide(F.b, F.c, P, F.a));
}

bool isInsideQuadrangle(Point &P, Quadrangle &F) {
    if (onOneSide(F.a, F.c, F.b, F.d))
        return (isInsideTriangle(P, Triangle(F.b, F.d, F.a)) || isInsideTriangle(P, Triangle(F.b, F.d, F.c)));
    else
        return (isInsideTriangle(P, Triangle(F.a, F.c, F.b)) || isInsideTriangle(P, Triangle(F.a, F.c, F.d)));
}

Intersection intersectLineWithPlane(Line L, Plane &G) {
    if (abs((L.b - L.a) % (G.b - G.a)) < EPS)
        return Intersection();
    long double t = -(G.d + (G.n * L.a)) / (G.n * (L.b - L.a));
    return Intersection(L.a + ((L.b - L.a) * t));
}

Intersection intersectRayWithPlane(Ray &R, Plane G) {
    Intersection I = intersectLineWithPlane(Line(R.a, R.b), G);
    if (I.hasIntersection && oneDirection(R.b - R.a, I.point - R.a))
        return Intersection(I.point);
    else
        return Intersection();
}

Intersection intersectRayWithTriangle(Ray &R, Triangle &F) {
    Intersection I = intersectRayWithPlane(R, Plane(F));
    if (I.hasIntersection && isInsideTriangle(I.point, F))
        return Intersection(I.point);
    else
        return Intersection();
}

Intersection intersectRayWithQuadrangle(Ray &R, Quadrangle &F) {
    if (onSameLine(F.a, F.b, F.c) || onSameLine(F.a, F.b, F.d) ||
        onSameLine(F.a, F.c, F.d) || onSameLine(F.b, F.c, F.d))
        return Intersection();
    Intersection I = intersectRayWithPlane(R, Plane(F.a, F.b, F.c));
    if (I.hasIntersection && isInsideQuadrangle(I.point, F))
        return Intersection(I.point);
    else
        return Intersection();
}

Intersection intersectRayWithSphere(Ray &R, Sphere &F) {
    Point P = projection(F.c, Line(R.a, R.b));
    if (dist(F.c, P) > F.r - EPS)
        return Intersection();
    Point X = P + (R.b - R.a) * sqrt(std::max(EPS, F.r * F.r - dist2(F.c, P))) / abs(R.b - R.a);
    Point Y = P - (R.b - R.a) * sqrt(std::max(EPS, F.r * F.r - dist2(F.c, P))) / abs(R.b - R.a);
    if (oneDirection(R.b - R.a, Y - R.a))
        return Intersection(Y);
    if (oneDirection(R.b - R.a, X - R.a))
        return Intersection(X);
    return Intersection();
}

Intersection Triangle::intersectWithRay(Ray &R) {
    return intersectRayWithTriangle(R, *this);
}

Intersection Quadrangle::intersectWithRay(Ray &R) {
    return intersectRayWithQuadrangle(R, *this);
}

Intersection Sphere::intersectWithRay(Ray &R) {
    return intersectRayWithSphere(R, *this);
}


Intersection intersectRayWithBBox(Point minBBox, Point maxBBox, Ray &R) {
    Point A(minBBox.x, minBBox.y, minBBox.z);
    Point B(minBBox.x, maxBBox.y, minBBox.z);
    Point C(maxBBox.x, maxBBox.y, minBBox.z);
    Point D(maxBBox.x, minBBox.y, minBBox.z);
    Point E(minBBox.x, minBBox.y, maxBBox.z);
    Point F(minBBox.x, maxBBox.y, maxBBox.z);
    Point G(maxBBox.x, maxBBox.y, maxBBox.z);
    Point H(maxBBox.x, minBBox.y, maxBBox.z);
    Quadrangle G1(A, B, C, D);
    Quadrangle G2(A, B, F, E);
    Quadrangle G3(B, C, G, F);
    Quadrangle G4(C, D, H, G);
    Quadrangle G5(D, A, E, H);
    Quadrangle G6(E, F, G, H);

    Intersection P1 = intersectRayWithQuadrangle(R, G1);
    Intersection P2 = intersectRayWithQuadrangle(R, G2);
    Intersection P3 = intersectRayWithQuadrangle(R, G3);
    Intersection P4 = intersectRayWithQuadrangle(R, G4);
    Intersection P5 = intersectRayWithQuadrangle(R, G5);
    Intersection P6 = intersectRayWithQuadrangle(R, G6);

    Point Closest = Point(INF, INF, INF);

    if (P1.hasIntersection && dist2(R.a, P1.point) < dist2(R.a, Closest))
        Closest = P1.point;
    if (P2.hasIntersection && dist2(R.a, P2.point) < dist2(R.a, Closest))
        Closest = P2.point;
    if (P3.hasIntersection && dist2(R.a, P3.point) < dist2(R.a, Closest))
        Closest = P3.point;
    if (P4.hasIntersection && dist2(R.a, P4.point) < dist2(R.a, Closest))
        Closest = P4.point;
    if (P5.hasIntersection && dist2(R.a, P5.point) < dist2(R.a, Closest))
        Closest = P5.point;
    if (P6.hasIntersection && dist2(R.a, P6.point) < dist2(R.a, Closest))
        Closest = P6.point;
    if (dist(Closest, Point(INF, INF, INF)) < EPS)
        return Intersection();
    return Intersection(Closest);

}

struct Light {
    Point p;
    long double e; //energy
    Light() {}
    Light(Point _p, long double _e) : p(_p), e(_e) {}
};


bool KDcompX(Figure *left, Figure *right) {
    return (left->getCenter().x < right->getCenter().x - EPS);
}

bool KDcompY(Figure *left, Figure *right) {
    return (left->getCenter().y < right->getCenter().y - EPS);
}

bool KDcompZ(Figure *left, Figure *right) {
    return (left->getCenter().z < right->getCenter().z - EPS);
}

struct KDNode {
    Figure *F;
    KDNode *left, *right;
    int n;
    Point minBBox, maxBBox;
    KDNode() : minBBox(Point(INF, INF, INF)), maxBBox(Point(-INF, -INF, -INF)), left(nullptr), right(nullptr), n(0) {}

    ~KDNode() {
        delete(left);
        delete(right);
    }

    void init(std::vector<Figure*>::iterator b, std::vector<Figure*>::iterator e, int coord) {
        n = e - b;
        //std::cout << n << "\n";

        if (!n)
            return;

        int m = n / 2;
        if (coord == 0)
            std::nth_element(b, b + m, e, KDcompX);
        else if (coord == 1)
            std::nth_element(b, b + m, e, KDcompY);
        else
            std::nth_element(b, b + m, e, KDcompZ);

        F = *(b + m);

        minBBox = F->getMinBBox();
        maxBBox = F->getMaxBBox();

        left = new KDNode();
        right = new KDNode();
        if (n > 1) {
            left->init(b, b + m, (coord + 1) % 3);
            right->init(b + m + 1, e, (coord + 1) % 3);

            minBBox.x = std::min(minBBox.x, std::min(left->minBBox.x, right->minBBox.x));
            minBBox.y = std::min(minBBox.y, std::min(left->minBBox.y, right->minBBox.y));
            minBBox.z = std::min(minBBox.z, std::min(left->minBBox.z, right->minBBox.z));
            maxBBox.x = std::max(maxBBox.x, std::max(left->maxBBox.x, right->maxBBox.x));
            maxBBox.y = std::max(maxBBox.y, std::max(left->maxBBox.y, right->maxBBox.y));
            maxBBox.z = std::max(maxBBox.z, std::max(left->maxBBox.z, right->maxBBox.z));
        }
    }

    Intersection intersectRayWithBBoxKDTree(Ray &R) {
        if (!n)
            return Intersection();
        return intersectRayWithBBox(minBBox, maxBBox, R);
    }

    std::pair<Point, Figure*> getClosestPointOnRay(Ray &R) { // гарантируем, что пересекается с BBox
        if (!n)
            return std::make_pair(Point(INF, INF, INF), nullptr);

        Intersection leftI = left->intersectRayWithBBoxKDTree(R);
        Intersection rightI = right->intersectRayWithBBoxKDTree(R);
        Intersection centerI = F->intersectWithRay(R);

        long double leftDist = leftI.hasIntersection ? dist(R.a, leftI.point) : INF;
        long double rightDist = rightI.hasIntersection ? dist(R.a, rightI.point) : INF;
        long double centerDist = centerI.hasIntersection ? dist(R.a, centerI.point) : INF;
        long double minDist = std::min(leftDist, rightDist);

        std::pair<Point, Figure*> closestLeft = std::make_pair(Point(INF, INF, INF), nullptr);
        std::pair<Point, Figure*> closestRight = std::make_pair(Point(INF, INF, INF), nullptr);
        std::pair<Point, Figure*> closestCenter = centerI.hasIntersection ? std::make_pair(centerI.point, F) :
                                                                            std::make_pair(Point(INF, INF, INF), nullptr);

        long double realDist = centerDist;
        if (abs(leftDist - minDist) < EPS) {
            if (realDist + EPS > leftDist && leftI.hasIntersection) {
                closestLeft = left->getClosestPointOnRay(R);
                leftDist = dist(R.a, closestLeft.first);
                if (realDist + EPS > leftDist)
                    realDist = leftDist;
            }
            if (realDist + EPS > rightDist && rightI.hasIntersection) {
                closestRight = right->getClosestPointOnRay(R);
                rightDist = dist(R.a, closestRight.first);
                if (realDist + EPS > rightDist)
                    realDist = rightDist;
            }
        } else if (abs(rightDist - minDist) < EPS) {
            if (realDist + EPS > rightDist && rightI.hasIntersection) {
                closestRight = right->getClosestPointOnRay(R);
                rightDist = dist(R.a, closestRight.first);
                if (realDist + EPS > rightDist)
                    realDist = rightDist;
            }
            if (realDist + EPS > leftDist && leftI.hasIntersection) {
                closestLeft = left->getClosestPointOnRay(R);
                leftDist = dist(R.a, closestLeft.first);
                if (realDist + EPS > leftDist)
                    realDist = leftDist;
            }
        }
/*
        if (leftI.hasIntersection)
            closestLeft = left->getClosestPointOnRay(R);
        if (rightI.hasIntersection)
            closestRight = right->getClosestPointOnRay(R);
        if (centerI.hasIntersection)
            closestCenter = std::make_pair(centerI.point, F);

        long double leftDist = dist(R.a, closestLeft.first);
        long double rightDist = dist(R.a, closestRight.first);
        long double centerDist = dist(R.a, closestCenter.first);

        realDist = std::min(std::min(leftDist, rightDist), centerDist);
*/
        if (abs(realDist - leftDist) < EPS)
            return closestLeft;
        if (abs(realDist - rightDist) < EPS)
            return closestRight;
        if (abs(realDist - centerDist) < EPS)
            return closestCenter;

        return std::make_pair(Point(INF, INF, INF), nullptr);
    }
};

struct KDTree {
    KDNode *root;

    KDTree() {}
    ~KDTree() {
        delete(root);
    }

    void init(std::vector<Figure*> &figures) {
        root = new KDNode();
        root->init(figures.begin(), figures.end(), 0);
    }

    std::pair<Point, Figure*> getClosestPointOnRay(Ray &R) {
        Intersection I = root->intersectRayWithBBoxKDTree(R);
        if (!I.hasIntersection)
            return std::make_pair(Point(INF, INF, INF), nullptr);
        return root->getClosestPointOnRay(R);
    }
};

struct GeomSpace {
    Point camera;
    Quadrangle *screen;
    std::vector<Figure*> figures;
    int width, height;
    KDTree kdTree;
    std::vector<Light> lights;

    GeomSpace() {}
    ~GeomSpace() {
        for (int i = 0; i < figures.size(); ++i)
            delete(figures[i]);
        delete screen;
    }

    void initKDTree() {
        kdTree.init(figures);
    }
};

#endif
