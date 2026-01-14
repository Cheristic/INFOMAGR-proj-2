#pragma once

#include "vec3.h"

struct photon
{
    point3 pos;
    char p[4];
    char phi, theta;
    short flag;
};

inline photon max(const photon& u, const photon& v)
{
    photon newPhoton;
    newPhoton.pos[0] = (u.pos[0] >= v.pos[0]) ? u.pos[0] : v.pos[0];
    newPhoton.pos[1] = (u.pos[1] >= v.pos[1]) ? u.pos[1] : v.pos[1];
    newPhoton.pos[2] = (u.pos[2] >= v.pos[2]) ? u.pos[2] : v.pos[2];
    return newPhoton;
}

inline photon min(const photon& u, const photon& v)
{
    photon newPhoton;
    newPhoton.pos[0] = (u.pos[0] <= v.pos[0]) ? u.pos[0] : v.pos[0];
    newPhoton.pos[1] = (u.pos[1] <= v.pos[1]) ? u.pos[1] : v.pos[1];
    newPhoton.pos[2] = (u.pos[2] <= v.pos[2]) ? u.pos[2] : v.pos[2];
    return newPhoton;
}