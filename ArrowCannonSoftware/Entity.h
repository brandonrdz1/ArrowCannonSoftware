#ifndef ENTITY_H
#define ENTITY_H

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <functional>
#include <stdexcept>
#include <cmath>

enum class direction { xDir, yDir, zDir };

class Entity {
public:
    double x, y, z, u, v, w;
    unsigned int amount = 1;
public:
    Entity();
    Entity(const std::string& inputString);
    Entity(double x, double y, double z);
    Entity(double x, double y, double z, double u, double v, double w);

    double getX() { return x; }
    double getY() { return y; }
    double getZ() { return z; }
    double getU() { return u; }
    double getV() { return v; }
    double getW() { return w; }

    void setX(double x) { this->x = x; }
    void setY(double y) { this->y = y; }
    void setZ(double z) { this->z = z; }
    void setU(double u) { this->u = u; }
    void setV(double v) { this->v = v; }
    void setW(double w) { this->w = w; }
    
    void print(bool velocity = false) {
        std::cout << "x: " << x << " y: " << y << " z: " << z;
        if (velocity) {
            std::cout << " | u: " << u << " v: " << v << " w: " << w;
        }
        std::cout << std::endl;
    }
    static std::string difference(Entity* pow, Entity* pro) {
        return ("dx: " + std::to_string(pro->x - pow->x) + " dy: " + std::to_string((pro->y + pro->getExplosionHeight()) - (pow->y+pow->getExplosionHeight())) + " dz: " + std::to_string(pro->z - pow->z));
    }
    double distance(Entity* a, Entity* b) {
        return sqrt((a->x - b->x) * (a->x - b->x) + (a->y - b->y) * (a->y - b->y) + (a->z - b->z) * (a->z - b->z));
    }
    virtual double* gainArray(Entity* powers, const int size, const double fraction = 1.0, const direction dir = direction::yDir) = 0;
    virtual void explosion(Entity* explosionSource, double exposure = 1.0) = 0;
    virtual void freefall(unsigned short int ticks) = 0;
    virtual double getExplosionHeight() const = 0;
    virtual ~Entity() = default;
};

class Tnt : public Entity {
public:
    static const double explosion_height;
    static const double gravity;
    static const double drag;

    Tnt();
    Tnt(const std::string& inputString);
    Tnt(double x, double y, double z);
    Tnt(double x, double y, double z, double u, double v, double w);

    double* gainArray(Entity* powers, const int size, const double fraction = 1.0, const direction dir = direction::yDir) override;
    void explosion(Entity* explosionSource, double exposure = 1.0) override;
    void freefall(unsigned short int ticks) override;
    double getExplosionHeight() const override;
};

class Arrow : public Entity {
public:
    static const double explosion_height;
    static const double gravity;
    static const double drag;

    Arrow();
    Arrow(const std::string& inputString);
    Arrow(double x, double y, double z);
    Arrow(double x, double y, double z, double u, double v, double w);

    double* gainArray(Entity* powers, const int size, const double fraction = 1.0, const direction dir = direction::yDir) override;
    void explosion(Entity* explosionSource, double exposure = 1.0) override;
    void freefall(unsigned short int ticks) override;
    double getExplosionHeight() const override;
};

#endif // ENTITY_H
