#ifndef HOMLIB_BACKEND_SEED_H_
#define HOMLIB_BACKEND_SEED_H_

// Seeds of the random generators of the GC-RANSAC and MAGSAC++ backends: those of GC-RANSAC's
// samplers, which MAGSAC++ uses too. GC-RANSAC seeds each of its generators from
// std::random_device; HomLib's copy of its uniform_random_generator.h (see CMakeLists.txt) calls
// next_generator() instead. Without a seed, next_generator() seeds from std::random_device too.

#include <cstdlib>
#include <optional>
#include <random>

namespace HomLib {
namespace backend_seed {

// The seed of the next generator created in this thread, if any
inline std::optional<unsigned int>& next_seed() {
    static thread_local std::optional<unsigned int> seed;
    return seed;
}

// A generator seeded with the next seed (the following generators get the next integers), or from
// std::random_device without a seed
inline std::mt19937 next_generator() {
    std::optional<unsigned int>& seed = next_seed();
    if (!seed) {
        std::random_device device;
        return std::mt19937(device());
    }
    return std::mt19937((*seed)++);
}

// In its scope, the generators created in this thread by next_generator() are seeded from seed
// (if any). With a seed, std::srand(seed) is called too, since GC-RANSAC also uses std::rand (in its
// adaptive reordering sampler, for one). std::rand is shared by the whole process, and std::srand
// can't be undone (the previous state can't be read): after a seeded run, std::rand continues from
// that seed. The same seed and settings give the same results on the same platform.
class Scope {
public:
    explicit Scope(std::optional<unsigned int> seed) : previous_(next_seed()) {
        next_seed() = seed;
        if (seed)
            std::srand(*seed);
    }
    ~Scope() { next_seed() = previous_; }
    Scope(const Scope&) = delete;
    Scope& operator=(const Scope&) = delete;

private:
    std::optional<unsigned int> previous_;
};

}  // namespace backend_seed
}  // namespace HomLib

#endif  // HOMLIB_BACKEND_SEED_H_
