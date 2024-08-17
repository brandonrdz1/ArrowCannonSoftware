#include "entity.h"

// Tnt's explosion height and drag
const double Tnt::explosion_height = 0.061250001192092896;
const double Tnt::gravity = -0.04;
const double Tnt::drag = 0.98;

// Arrow's explosion height and drag
const double Arrow::explosion_height = 0.12999999523162842;
const double Arrow::gravity = -0.05f;
const double Arrow::drag = 0.99f;

Entity::Entity() : x(0), y(0), z(0), u(0), v(0), w(0) {}

Entity::Entity(const std::string& inputString) {
    std::istringstream(inputString) >> x >> y >> z >> u >> v >> w;
}

Entity::Entity(double x, double y, double z) : x(x), y(y), z(z), u(0), v(0), w(0) {}

Entity::Entity(double x, double y, double z, double u, double v, double w) : x(x), y(y), z(z), u(u), v(v), w(w) {}

// Tnt subclass implementations
Tnt::Tnt() : Entity() {}

Tnt::Tnt(const std::string& inputString) : Entity(inputString) {}

Tnt::Tnt(double x, double y, double z) : Entity(x, y, z) {}

Tnt::Tnt(double x, double y, double z, double u, double v, double w) : Entity(x, y, z, u, v, w) {}

double* Tnt::gainArray(Entity* powers, const int size, const double fraction, const direction dir) {
    if (size == 0) throw std::runtime_error("No tnts added!");

    double* gainArr = new double[size];
    double dx, dy, dz, dist, f;

    for (int i = 0; i < size; ++i) {
        Entity& booster = powers[i];
        dx = this->x - booster.x;
        dy = this->y - (booster.y + Tnt::explosion_height);
        dz = this->z - booster.z;
        dist = std::sqrt(dx * dx + dy * dy + dz * dz);

        if (dist > 8.0 || dist == 0.0) {
            gainArr[i] = 0.0;
            continue;
        }

        f = (1.0 - dist / 8.0) * fraction / dist;

        switch (dir) {
        case direction::xDir:
            gainArr[i] = dx * f;
            break;
        case direction::yDir:
            gainArr[i] = dy * f;
            break;
        case direction::zDir:
            gainArr[i] = dz * f;
            break;
        }
    }
    return gainArr;
}

void Tnt::explosion(Entity* explosionSource, double exposure) {
    // Implement the explosion logic here
    double dx, dy, dz, dist, f;
    dx = this->x - explosionSource->x;
    dy = this->y - (explosionSource->y + Tnt::explosion_height);
    dz = this->z - explosionSource->z;

    dist = std::sqrt(dx * dx + dy * dy + dz * dz);
    if (dist > 8.0 || dist == 0.0) return;
    f = (1.0 - dist / 8.0) * explosionSource->amount / dist * exposure;

    this->u += dx * f;
    this->v += dy * f;
    this->w += dz * f;
}

void Tnt::freefall(unsigned short int ticks) {
    // Update position based on current velocity and apply drag
    for (unsigned short int i = 0; i < ticks; ++i) {
        this->v += Tnt::gravity;

        this->x += this->u;
        this->y += this->v;
        this->z += this->w;

        // Apply drag
        this->u *= Tnt::drag;
        this->v *= Tnt::drag;
        this->w *= Tnt::drag;
    }
}

double Tnt::getExplosionHeight() const {
    return Tnt::explosion_height;
}

// Arrow subclass implementations
Arrow::Arrow() : Entity() {}

Arrow::Arrow(const std::string& inputString) : Entity(inputString) {}

Arrow::Arrow(double x, double y, double z) : Entity(x, y, z) {}

Arrow::Arrow(double x, double y, double z, double u, double v, double w) : Entity(x, y, z, u, v, w) {}

double* Arrow::gainArray(Entity* powers, const int size, const double fraction, const direction dir) {
    if (size == 0) throw std::runtime_error("No arrows added!");

    double* gainArr = new double[size];
    double dx, dy, dz, dist, f, directionMagnitude;

    for (int i = 0; i < size; ++i) {
        Entity& booster = powers[i];
        dx = this->x - booster.x;
        dy = this->y - (booster.y + this->getExplosionHeight()); // FXME
        dz = this->z - booster.z;
        dist = std::sqrt(dx * dx + dy * dy + dz * dz);
        std::cout << "dist: " << dist << std::endl;
        std::cout << "dy: " << dy << std::endl;

        if (dist > 8.0 || dist == 0.0) {
            gainArr[i] = 0.0;
            continue;
        }

        f = (1.0 - dist / 8.0) * fraction;
        directionMagnitude = std::sqrt(dx * dx + dy * dy + dz * dz);

        switch (dir) {
        case direction::xDir:
            gainArr[i] = dx / directionMagnitude * f;
            break;
        case direction::yDir:
            gainArr[i] = dy / directionMagnitude * f;
            break;
        case direction::zDir:
            gainArr[i] = dz / directionMagnitude * f;
            break;
        }
    }

    return gainArr;
}

void Arrow::explosion(Entity* explosionSource, double exposure) {
    double dx, dy, dz, dist, f, dirMag;
    dx = this->x - explosionSource->x;
    dy = this->y - (explosionSource->y + Tnt::explosion_height);
    dz = this->z - explosionSource->z;

    dist = std::sqrt(dx * dx + dy * dy + dz * dz);
    if (dist > 8.0 || dist == 0.0) return;
    f = (1.0 - dist / 8.0) * explosionSource->amount * exposure;
    dy = (this->y + Arrow::explosion_height) - (explosionSource->y + Tnt::explosion_height); // Eye Height
    dirMag = std::sqrt(dx * dx + dy * dy + dz * dz);

    this->u += dx / dirMag * f;
    this->v += dy / dirMag * f;
    this->w += dz / dirMag * f;
}

void Arrow::freefall(unsigned short int ticks) {
    // Update position based on current velocity and apply drag
    for (unsigned short int i = 0; i < ticks; ++i) {
        this->x += this->u;
        this->y += this->v;
        this->z += this->w;

        // Apply drag
        this->u *= Arrow::drag;
        this->v *= Arrow::drag;
        this->w *= Arrow::drag;

        this->v += Arrow::gravity;
    }
}
/*
vector<double> powerLocations;
    double gravity = -0.05f;
    double drag = 0.99f;
    cout << "start: " << arrow.pos[1] << "\t" << arrow.mot[1] << endl;
    for (int i = 0; i < ticks; i++) {
        arrow.pos[1] += arrow.mot[1];
        arrow.mot[1] *= drag;
        arrow.mot[1] += gravity;
        cout << i << " " << arrow.pos[1] << "\t" << arrow.mot[1] << endl;
        powerLocations.push_back(arrow.pos[1] + 0.12999999523162842 - Entity::explosion_height);
    }

    for (int i = 0; i < ticks; i++) {
        cout << i << " " << powerLocations[i] << endl;
    }*/
double Arrow::getExplosionHeight() const {
    return Arrow::explosion_height;
}
