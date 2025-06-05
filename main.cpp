#include <iostream>
#include <fstream>
#include <vector>
#include <glm/glm.hpp>

// Structure for material properties (Phong model)
struct Material {
    glm::vec3 ambient;
    glm::vec3 diffuse;
    glm::vec3 specular;
    float shininess;
    float reflectivity;
};

// Structure for light properties
struct Light {
    glm::vec3 position;
    glm::vec3 ambient;
    glm::vec3 diffuse;
    glm::vec3 specular;
};

// Structure for a sphere
struct Sphere {
    glm::vec3 center;
    float radius;
    Material material;
};

// Structure for a plane
struct Plane {
    glm::vec3 point;  // A point on the plane (e.g., at y=0)
    glm::vec3 normal; // Plane normal (e.g., (0,1,0) for XZ plane)
    Material material;
};

// Ray structure
struct Ray {
    glm::vec3 origin;
    glm::vec3 direction;
};

// Intersection result
struct Intersection {
    bool hit;
    glm::vec3 point;
    glm::vec3 normal;
    float t;
    Material material; // Store material of intersected object
};

// Compute ray-sphere intersection
Intersection intersectSphere(const Ray& ray, const Sphere& sphere) {
    Intersection isect = {false, glm::vec3(0), glm::vec3(0), 0.0f, sphere.material};
    glm::vec3 oc = ray.origin - sphere.center;
    float a = glm::dot(ray.direction, ray.direction);
    float b = 2.0f * glm::dot(oc, ray.direction);
    float c = glm::dot(oc, oc) - sphere.radius * sphere.radius;
    float discriminant = b * b - 4 * a * c;

    if (discriminant < 0) return isect;

    float t = (-b - glm::sqrt(discriminant)) / (2.0f * a);
    if (t < 0) {
        t = (-b + glm::sqrt(discriminant)) / (2.0f * a);
        if (t < 0) return isect;
    }

    isect.hit = true;
    isect.t = t;
    isect.point = ray.origin + t * ray.direction;
    isect.normal = glm::normalize(isect.point - sphere.center);
    return isect;
}

// Compute ray-plane intersection
Intersection intersectPlane(const Ray& ray, const Plane& plane) {
    Intersection isect = {false, glm::vec3(0), plane.normal, 0.0f, plane.material};
    float denom = glm::dot(plane.normal, ray.direction);
    if (glm::abs(denom) > 1e-6) { // Avoid division by zero
        glm::vec3 p0 = plane.point - ray.origin;
        float t = glm::dot(p0, plane.normal) / denom;
        if (t >= 0) {
            isect.hit = true;
            isect.t = t;
            isect.point = ray.origin + t * ray.direction;
            return isect;
        }
    }
    return isect;
}

// Compute Phong shading
glm::vec3 computePhong(const Intersection& isect, const Light& light, const glm::vec3& viewPos) {
    // Ambient
    glm::vec3 ambient = light.ambient * isect.material.ambient;

    // Diffuse
    glm::vec3 norm = glm::normalize(isect.normal);
    glm::vec3 lightDir = glm::normalize(light.position - isect.point);
    float diff = glm::max(glm::dot(norm, lightDir), 0.0f);
    glm::vec3 diffuse = light.diffuse * (diff * isect.material.diffuse);

    // Specular
    glm::vec3 viewDir = glm::normalize(viewPos - isect.point);
    glm::vec3 reflectDir = glm::reflect(-lightDir, norm);
    float spec = glm::pow(glm::max(glm::dot(viewDir, reflectDir), 0.0f), isect.material.shininess);
    glm::vec3 specular = light.specular * (spec * isect.material.specular);

    return ambient + diffuse + specular;
}

glm::vec3 randomInHemisphereUniform(const glm::vec3& normal) {
    float xi1 = static_cast<float>(rand()) / RAND_MAX;
    float xi2 = static_cast<float>(rand()) / RAND_MAX;

    float r = sqrt(1.0f - xi1 * xi1);
    float phi = 2.0f * M_PI * xi2;
    float x = r * cos(phi);
    float y = r * sin(phi);
    float z = xi1;

    // Convert from tangent space to world space
    glm::vec3 up = fabs(normal.z) < 0.999f ? glm::vec3(0, 0, 1) : glm::vec3(1, 0, 0);
    glm::vec3 tangent = glm::normalize(glm::cross(up, normal));
    glm::vec3 bitangent = glm::cross(normal, tangent);

    return glm::normalize(x * tangent + y * bitangent + z * normal);
}

glm::vec3 traceRay(const Ray& ray, const Sphere& sphere1, const Sphere& sphere2,
                   const std::vector<Plane>& roomWalls, const Light& light,
                   const glm::vec3& viewPos, int depth = 0) {
    const int NUM_SAMPLES = 16;

    if (depth > 3) return glm::vec3(0.0f); // Max recursion depth

    Intersection sphere1Isect = intersectSphere(ray, sphere1);
    Intersection sphere2Isect = intersectSphere(ray, sphere2);
    Intersection closest = {false, glm::vec3(0), glm::vec3(0), std::numeric_limits<float>::max(), Material{}};

    // Check sphere intersections
    if (sphere1Isect.hit && sphere1Isect.t < closest.t) closest = sphere1Isect;
    if (sphere2Isect.hit && sphere2Isect.t < closest.t) closest = sphere2Isect;

    // Check room wall intersections
    for (const auto& wall : roomWalls) {
        Intersection wallIsect = intersectPlane(ray, wall);
        if (wallIsect.hit && wallIsect.t < closest.t) closest = wallIsect;
    }

    if (!closest.hit) return glm::vec3(0.1f); // Background

    glm::vec3 color = computePhong(closest, light, viewPos);

    // Reflections
    if (closest.material.reflectivity > 0.0f) {
        glm::vec3 reflectDir = glm::reflect(ray.direction, closest.normal);
        Ray reflectRay = {closest.point + 0.001f * closest.normal, reflectDir};
        glm::vec3 reflected = traceRay(reflectRay, sphere1, sphere2, roomWalls, light, viewPos, depth + 1);
        color = glm::mix(color, reflected, closest.material.reflectivity);
    }

    // Monte Carlo diffuse global illumination
    if (depth < 2 && closest.material.reflectivity < 1.0f) {
        glm::vec3 indirectLighting(0.0f);
        for (int i = 0; i < NUM_SAMPLES; ++i) {
            glm::vec3 bounceDir = randomInHemisphereUniform(closest.normal);
            Ray bounceRay = {closest.point + 0.001f * closest.normal, bounceDir};
            glm::vec3 bounceColor = traceRay(bounceRay, sphere1, sphere2, roomWalls, light, viewPos, depth + 1);
            indirectLighting += bounceColor;
        }
        indirectLighting /= static_cast<float>(NUM_SAMPLES);
        color += 0.5f * indirectLighting * closest.material.diffuse;
    }

    return glm::clamp(color, 0.0f, 1.0f);
}

int main() {
    // Image settings
    const int width = 1920;
    const int height = 1080;
    std::vector<glm::vec3> framebuffer(width * height);

    // Camera settings
    glm::vec3 cameraPos(0.0f, 1.0f, 6.0f);
    float fov = glm::radians(45.0f);
    float aspect = static_cast<float>(width) / height;
    float tanFov = glm::tan(fov / 2.0f);

    // Scene setup: first sphere
    Sphere sphere1;
    sphere1.center = glm::vec3(-1.0f, 1.0f, 0.0f);
    sphere1.radius = 1.0f;
    sphere1.material = {
        glm::vec3(1.0f, 0.5f, 0.31f), // ambient
        glm::vec3(1.0f, 0.5f, 0.31f), // diffuse
        glm::vec3(1.0f, 1.0f, 1.0f),  // specular
        128.0f,                       // shininess
        0.4f                          // more reflective
    };

    // Scene setup: second sphere
    Sphere sphere2;
    sphere2.center = glm::vec3(2.0f, 1.0f, 0.0f); // Positioned to the right
    sphere2.radius = 1.0f;
    sphere2.material = {
        glm::vec3(0.3f, 0.8f, 0.3f),  // ambient (greenish)
        glm::vec3(0.3f, 0.8f, 0.3f),  // diffuse
        glm::vec3(1.0f, 1.0f, 1.0f),  // specular
        128.0f,                       // shininess
        0.4f                          // more reflective
    };

    // Room material
    Material roomMaterial = {
        glm::vec3(0.7f, 0.7f, 0.8f),  // ambient (light blue-grey)
        glm::vec3(0.7f, 0.7f, 0.8f),  // diffuse
        glm::vec3(1.0f, 1.0f, 1.0f),  // specular
        64.0f,                        // shininess
        0.15f                         // less reflective
    };

    // Create room walls
    std::vector<Plane> roomWalls;

    // Floor (y = 0)
    Plane floor;
    floor.point = glm::vec3(0.0f, 0.0f, 0.0f);
    floor.normal = glm::vec3(0.0f, 1.0f, 0.0f);
    floor.material = roomMaterial;
    floor.material.reflectivity = 0.35f; // Less reflective floor
    roomWalls.push_back(floor);

    // Ceiling (y = 4)
    Plane ceiling;
    ceiling.point = glm::vec3(0.0f, 4.0f, 0.0f);
    ceiling.normal = glm::vec3(0.0f, -1.0f, 0.0f);
    ceiling.material = roomMaterial;
    roomWalls.push_back(ceiling);

    // Left wall (x = -4)
    Plane leftWall;
    leftWall.point = glm::vec3(-4.0f, 0.0f, 0.0f);
    leftWall.normal = glm::vec3(1.0f, 0.0f, 0.0f);
    leftWall.material = roomMaterial;
    leftWall.material.ambient = glm::vec3(0.8f, 0.6f, 0.6f); // Slight red tint
    leftWall.material.diffuse = glm::vec3(0.8f, 0.6f, 0.6f);
    roomWalls.push_back(leftWall);

    // Right wall (x = 4)
    Plane rightWall;
    rightWall.point = glm::vec3(4.0f, 0.0f, 0.0f);
    rightWall.normal = glm::vec3(-1.0f, 0.0f, 0.0f);
    rightWall.material = roomMaterial;
    rightWall.material.ambient = glm::vec3(0.6f, 0.8f, 0.6f); // Slight green tint
    rightWall.material.diffuse = glm::vec3(0.6f, 0.8f, 0.6f);
    roomWalls.push_back(rightWall);

    // Back wall (z = -4)
    Plane backWall;
    backWall.point = glm::vec3(0.0f, 0.0f, -4.0f);
    backWall.normal = glm::vec3(0.0f, 0.0f, 1.0f);
    backWall.material = roomMaterial;
    roomWalls.push_back(backWall);

    // Front wall (z = 8) - behind camera, won't be visible but completes the room
    Plane frontWall;
    frontWall.point = glm::vec3(0.0f, 0.0f, 8.0f);
    frontWall.normal = glm::vec3(0.0f, 0.0f, -1.0f);
    frontWall.material = roomMaterial;
    roomWalls.push_back(frontWall);

    // Light setup - dimmer light positioned inside the room
    Light light;
    light.position = glm::vec3(1.2f, 3.0f, 2.0f);
    light.ambient = glm::vec3(0.15f);
    light.diffuse = glm::vec3(0.6f);
    light.specular = glm::vec3(0.8f);

    // Ray tracing loop
    for (int j = 0; j < height; ++j) {
        for (int i = 0; i < width; ++i) {
            // Compute normalized device coordinates
            float x = (2.0f * i / width - 1.0f) * aspect * tanFov;
            float y = (1.0f - 2.0f * j / height) * tanFov;
            glm::vec3 dir = glm::normalize(glm::vec3(x, y, -1.0f));

            // Create ray
            Ray ray = {cameraPos, dir};

            // Trace ray through the scene
            glm::vec3 color = traceRay(ray, sphere1, sphere2, roomWalls, light, cameraPos);

            // Store in framebuffer
            framebuffer[j * width + i] = color;
        }
    }

    // Output to PPM file
    std::ofstream ofs("output.ppm");
    ofs << "P3\n" << width << " " << height << "\n255\n";
    for (const auto& color : framebuffer) {
        int r = static_cast<int>(color.r * 255.99f);
        int g = static_cast<int>(color.g * 255.99f);
        int b = static_cast<int>(color.b * 255.99f);
        ofs << r << " " << g << " " << b << "\n";
    }
    ofs.close();

    std::cout << "Ray tracing complete. Output saved to output.ppm\n";
    return 0;
}