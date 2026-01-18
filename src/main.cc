//==============================================================================================
// Originally written in 2016 by Peter Shirley <ptrshrl@gmail.com>
//
// To the extent possible under law, the author(s) have dedicated all copyright and related and
// neighboring rights to this software to the public domain worldwide. This software is
// distributed without any warranty.
//
// You should have received a copy (see file COPYING.txt) of the CC0 Public Domain Dedication
// along with this software. If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.
//==============================================================================================

#include "rtweekend.h"

#include "bvh.h"
#include "camera.h"
#include "hittable.h"
#include "hittable_list.h"
#include "material.h"
#include "quad.h"
#include "sphere.h"
#include "texture.h"
#include "light.h"
#include "photon_mapping/kdtree.h"

void cornell_box() {
    hittable_list world;


    auto red   = make_shared<lambertian>(color(.65, .05, .05));
    auto white = make_shared<lambertian>(color(.73, .73, .73));
    auto green = make_shared<lambertian>(color(.12, .45, .15));
    auto light_mat = make_shared<diffuse_light>(color(30, 30, 30));

    world.add(make_shared<quad>(point3(555,0,0), vec3(0,555,0), vec3(0,0,555), green));
    world.add(make_shared<quad>(point3(0,0,0), vec3(0,555,0), vec3(0,0,555), red));
    const std::shared_ptr<light> light_quad = make_shared<light>(point3(378, 554, 332), vec3(-200, 0, 0), vec3(0, 0, -105), light_mat, vec3(1,1,1));
    world.add(light_quad);
    world.add(make_shared<quad>(point3(0,0,0), vec3(555,0,0), vec3(0,0,555), white));
    world.add(make_shared<quad>(point3(555,555,555), vec3(-555,0,0), vec3(0,0,-555), white));
    world.add(make_shared<quad>(point3(0,0,555), vec3(555,0,0), vec3(0,555,0), white));

    shared_ptr<hittable> box1 = box(point3(0,0,0), point3(200,200,200), white);
    box1 = make_shared<rotate_y>(box1, 55);
    box1 = make_shared<translate>(box1, vec3(265,0,180));
    world.add(box1);

    world.add(make_shared<sphere>(
        point3(120, 90, 280), 90, make_shared<metal>(color(0.8, 0.8, 0.9), 0.8)));
    camera cam;

    cam.aspect_ratio      = 16.0 / 9.0;
    cam.image_width       = 1200;
    cam.samples_per_pixel = 4;
    cam.max_depth         = 100;
    cam.background        = color(0,0,0);

    cam.vfov     = 40;
    cam.lookfrom = point3(278, 278, -800);
    cam.lookat   = point3(278, 278, 0);
    cam.vup      = vec3(0,1,0);

    cam.defocus_angle = 0;

    PhotonMap map;

    map.nPhotonsGlobal = 20000; // should be 100000
    map.maxDepth = 100;
    map.nPhotonsCaustic = map.nPhotonsGlobal * 100;
    map.nEstimationPhotons = 100;
    map.finalGatheringDepth = 0;

    map.build(world, light_quad);

    //cam.render_photons(world, map);
    cam.render(world, map);
}

void new_cornell_box() {
    hittable_list world;


    auto red = make_shared<lambertian>(color(.65, .05, .05));
    auto white = make_shared<lambertian>(color(.73, .73, .73));
    auto green = make_shared<lambertian>(color(.12, .45, .15));
    auto light_mat = make_shared<diffuse_light>(color(30, 30, 30));

    world.add(make_shared<quad>(point3(555, 0, 0), vec3(0, 555, 0), vec3(0, 0, 555), green));
    world.add(make_shared<quad>(point3(0, 0, 0), vec3(0, 555, 0), vec3(0, 0, 555), red));
    const std::shared_ptr<light> light_quad = make_shared<light>(point3(378, 554, 332), vec3(-200, 0, 0), vec3(0, 0, -105), light_mat, vec3(1, 1, 1));
    world.add(light_quad);
    world.add(make_shared<quad>(point3(0, 0, 0), vec3(555, 0, 0), vec3(0, 0, 555), white));
    world.add(make_shared<quad>(point3(555, 555, 555), vec3(-555, 0, 0), vec3(0, 0, -555), white));
    world.add(make_shared<quad>(point3(0, 0, 555), vec3(555, 0, 0), vec3(0, 555, 0), white));

    // back wall
    world.add(make_shared<quad>(point3(0, 0, -20), vec3(555, 0, 0), vec3(0, 555, 0), white));

    shared_ptr<hittable> box1 = box(point3(0, 0, 0), point3(200, 200, 200), white);
    box1 = make_shared<rotate_y>(box1, 55);
    box1 = make_shared<translate>(box1, vec3(265, 0, 180));
    world.add(box1);

    world.add(make_shared<sphere>(
        point3(120, 90, 280), 90, make_shared<metal>(color(0.8, 0.8, 0.9), 0.8)));
    camera cam;

    cam.aspect_ratio = 1;
    cam.image_width = 1200;
    cam.samples_per_pixel = 4;
    cam.max_depth = 100;
    cam.background = color(0, 0, 0);

    cam.vfov = 100;
    cam.lookfrom = point3(278, 278, -1);
    cam.lookat = point3(278, 278, 0);
    cam.vup = vec3(0, 1, 0);

    cam.defocus_angle = 0;

    PhotonMap map;

    map.nPhotonsGlobal = 20000; // should be 100000
    map.maxDepth = 100;
    map.nPhotonsCaustic = map.nPhotonsGlobal * 100;
    map.nEstimationPhotons = 100;
    map.finalGatheringDepth = 0;

    map.build(world, light_quad);

    //cam.render_photons(world, map);
    cam.render(world, map);
}

void empty()
{
    hittable_list world;

    auto red = make_shared<lambertian>(color(.65, .05, .05));
    auto white = make_shared<lambertian>(color(.73, .73, .73));
    auto green = make_shared<lambertian>(color(.12, .45, .15));
    auto light_mat = make_shared<diffuse_light>(color(30, 30, 30));

    world.add(make_shared<quad>(point3(555, 0, 0), vec3(0, 555, 0), vec3(0, 0, 555), green));
    world.add(make_shared<quad>(point3(0, 0, 0), vec3(0, 555, 0), vec3(0, 0, 555), red));
    const std::shared_ptr<light> light_quad = make_shared<light>(point3(378, 554, 332), vec3(-200, 0, 0), vec3(0, 0, -105), light_mat, vec3(1, 1, 1));
    world.add(light_quad);
    world.add(make_shared<quad>(point3(0, 0, 0), vec3(555, 0, 0), vec3(0, 0, 555), white));
    world.add(make_shared<quad>(point3(555, 555, 555), vec3(-555, 0, 0), vec3(0, 0, -555), white));
    world.add(make_shared<quad>(point3(0, 0, 555), vec3(555, 0, 0), vec3(0, 555, 0), white));

    camera cam;

    cam.aspect_ratio = 16.0 / 9.0;
    cam.image_width = 1200;
    cam.samples_per_pixel = 100;
    cam.max_depth = 100;
    cam.background = color(0, 0, 0);

    cam.vfov = 40;
    cam.lookfrom = point3(278, 278, -800);
    cam.lookat = point3(278, 278, 0);
    cam.vup = vec3(0, 1, 0);

    cam.defocus_angle = 0;

    PhotonMap map;

    map.nPhotonsGlobal = 1000000;
    map.maxDepth = 100;
    map.nPhotonsCaustic = map.nPhotonsGlobal * 100;
    map.nEstimationPhotons = 1;
    map.finalGatheringDepth = 0;

    map.build(world, light_quad);

    cam.render(world, map);
}

void visualize_photon_map() {
    hittable_list world;


    auto red = make_shared<lambertian>(color(.65, .05, .05));
    auto white = make_shared<lambertian>(color(.73, .73, .73));
    auto green = make_shared<lambertian>(color(.12, .45, .15));
    auto light_mat = make_shared<diffuse_light>(color(30, 30, 30));

    world.add(make_shared<quad>(point3(555, 0, 0), vec3(0, 555, 0), vec3(0, 0, 555), green));
    world.add(make_shared<quad>(point3(0, 0, 0), vec3(0, 555, 0), vec3(0, 0, 555), red));
    const std::shared_ptr<light> light_quad = make_shared<light>(point3(378, 554, 332), vec3(-200, 0, 0), vec3(0, 0, -105), light_mat, vec3(1, 1, 1));
    world.add(light_quad);
    world.add(make_shared<quad>(point3(0, 0, 0), vec3(555, 0, 0), vec3(0, 0, 555), white));
    world.add(make_shared<quad>(point3(555, 555, 555), vec3(-555, 0, 0), vec3(0, 0, -555), white));
    world.add(make_shared<quad>(point3(0, 0, 555), vec3(555, 0, 0), vec3(0, 555, 0), white));

    camera cam;

    cam.aspect_ratio = 16.0 / 9.0;
    cam.image_width = 600;
    cam.samples_per_pixel = 1;
    cam.max_depth = 100;
    cam.background = color(0, 0, 0);

    cam.vfov = 40;
    cam.lookfrom = point3(278, 278, -800);
    cam.lookat = point3(278, 278, 0);
    cam.vup = vec3(0, 1, 0);

    cam.defocus_angle = 0;

    PhotonMap map;

    map.nPhotonsGlobal = 100000;
    map.maxDepth = 100;
    map.nPhotonsCaustic = map.nPhotonsGlobal * 100;
    map.nEstimationPhotons = 100;
    map.finalGatheringDepth = 0;

    //std::clog << "beginning build";
    const clock_t begin_time = clock();
    // do something

    map.build(world, light_quad);
    //std::clog << "finished build in " << float(clock() - begin_time) / CLOCKS_PER_SEC;

    cam.render_photons(world, map);
}
int main() {
    switch (1)
    {
    case 2:
        visualize_photon_map();
        break;
    case 1:
        cornell_box();
        break;
    case 3:
        new_cornell_box();
        break;
    default:
        empty();
        break;
    }
    
}
