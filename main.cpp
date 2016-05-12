#ifndef _MAIN
#define _MAIN

#include "geomspace.h"
#include "draw.h"
#include <fstream>
#include <iostream>
#include <string>
#include <iomanip>

GeomSpace buildInputSpace() {
    GeomSpace gs;

    std::ifstream in;
    in.open("geomSpace.in", std::ios::in);

    std::string object;
    while (in >> object) {
        if (object == "end") {
            break;
        } else if (object == "size") {
            in >> gs.width >> gs.height;
        } else if (object == "camera") {
            in >> gs.camera;
        } else if (object == "screen") {
            Point a;
            long double w, h;
            in >> a >> w >> h;
            gs.screen = new Quadrangle(a, w, h, gs.camera);
            //std::cout << gs.screen->a << "\n"<< gs.screen->b << "\n"<< gs.screen->c << "\n"<< gs.screen->d << "\n";
        } else if (object == "Triangle") {
            Point a, b, c;
            in >> a >> b >> c;
            long double alpha = 0.0;
            in >> alpha;
            std::string typeColor;
            in >> typeColor;
            if (typeColor == "color") {
                Color clr;
                in >> clr;
                gs.figures.push_back(new Triangle(a, b, c, alpha, clr));
            } else {
                Texture* texture = new Texture();
                in >> *texture;
                gs.figures.push_back(new Triangle(a, b, c, alpha, texture));
            }
        } else if (object == "Quadrangle") {
            Point a, b, c, d;
            in >> a >> b >> c >> d;
            long double alpha = 0.0;
            in >> alpha;
            std::string typeColor;
            in >> typeColor;
            if (typeColor == "color") {
                Color clr;
                in >> clr;
                gs.figures.push_back(new Quadrangle(a, b, c, d, alpha, clr));
            } else {
                Texture* texture = new Texture();
                in >> *texture;
                gs.figures.push_back(new Quadrangle(a, b, c, d, alpha, texture));
            }
        } else if (object == "Sphere") {
            Point c;
            long double r;
            in >> c >> r;
            long double alpha = 0.0;
            in >> alpha;
            std::string typeColor;
            in >> typeColor;
            Color clr;
            in >> clr;
            gs.figures.push_back(new Sphere(c, r, alpha, clr));
        } else if (object == "Light") {
            Point c;
            long double i;
            in >> c >> i;
            gs.lights.push_back(Light(c, i));
        } else if (object == "File") {
            std::string fileName;
            in >> fileName;
            long double alpha = 0.0;
            in >> alpha;
            std::string typeColor;
            in >> typeColor;
            Color clr;
            Texture* texture = new Texture();
            if (typeColor == "color") {
                in >> clr;
            } else {
                in >> *texture;
            }
            std::ifstream fin;
            fin.open(fileName, std::ios::in);
            std::string s;
            while (fin >> s) {
                if (s == "loop") {
                    Point a, b, c;
                    std::string p;
                    fin >> p >> a;
                    fin >> p >> b;
                    fin >> p >> c;
                    if (typeColor == "color")
                        gs.figures.push_back(new Triangle(a, b, c, alpha, clr));
                    else
                        gs.figures.push_back(new Triangle(a, b, c, alpha, texture));
                } else if (s == "endsolid") {
                    break;
                }
            }
        }
    }
    gs.initKDTree();
    in.close();

    return gs;
}

int main() {
    GeomSpace geomSpace = buildInputSpace();
    drawThePicture(geomSpace);
    return 0;
}


#endif
