#ifndef NOISE_PLANE_FIT_HPP
#define NOISE_PLANE_FIT_HPP

#include "LensingConfig.hpp"
#include "UniversalUtils.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

namespace NoisePlaneFit {

    // ==========================================
    // Function: Fit the Stage-3 outer-noise background plane
    // Method: Sample the in-chip square shell between the configured outer and inner sizes,
    //         retaining only finite weight-one pixels from the initial source amplifier.
    // ==========================================
    inline bool fitNoiseRegionPlane(
        const std::vector<float>& localImage,
        const std::vector<int>& localWeight,
        int regionSize,
        int innerSize,
        int sourceOffset,
        int localStartX,
        int localStartY,
        int chipWidth,
        int chipHeight,
        int ccdSplit,
        int sourceChipX,
        double& aa,
        double& bb,
        double& cc) {
        if (regionSize <= 0 || innerSize <= 0 || innerSize >= regionSize
            || (regionSize - innerSize) % 2 != 0 || sourceOffset < 0
            || chipWidth <= 0 || chipHeight <= 0
            || (ccdSplit != 1 && ccdSplit != 2)
            || sourceChipX < 0 || sourceChipX >= chipWidth) {
            return false;
        }

        const std::size_t regionSide = static_cast<std::size_t>(regionSize);
        const std::size_t expectedSize = regionSide * regionSide;
        if (localImage.size() != expectedSize || localWeight.size() != expectedSize) {
            return false;
        }

        const int innerStart = (regionSize - innerSize) / 2;
        const int innerEnd = innerStart + innerSize;
        const int amplifierBoundary = chipWidth / 2;
        const int sourceAmplifier = sourceChipX < amplifierBoundary ? 0 : 1;

        std::vector<Point3D> points;
        points.reserve(expectedSize
                       - static_cast<std::size_t>(innerSize)
                             * static_cast<std::size_t>(innerSize));
        std::size_t geometricCandidates = 0U;

        for (int localY = 0; localY < regionSize; ++localY) {
            for (int localX = 0; localX < regionSize; ++localX) {
                const bool insideInner = localX >= innerStart && localX < innerEnd
                                      && localY >= innerStart && localY < innerEnd;
                if (insideInner) {
                    continue;
                }

                const int chipX = localStartX + localX;
                const int chipY = localStartY + localY;
                if (chipX < 0 || chipX >= chipWidth || chipY < 0 || chipY >= chipHeight) {
                    continue;
                }
                if (ccdSplit == 2) {
                    const int pixelAmplifier = chipX < amplifierBoundary ? 0 : 1;
                    if (pixelAmplifier != sourceAmplifier) {
                        continue;
                    }
                }

                ++geometricCandidates;
                const std::size_t localIndex = static_cast<std::size_t>(localY) * regionSide
                                             + static_cast<std::size_t>(localX);
                if (localWeight[localIndex] != 1
                    || !std::isfinite(localImage[localIndex])) {
                    continue;
                }

                points.push_back({
                    static_cast<double>(localX - sourceOffset + 1),
                    static_cast<double>(localY - sourceOffset + 1),
                    static_cast<double>(localImage[localIndex])
                });
            }
        }

        if (geometricCandidates == 0U) {
            return false;
        }
        const std::size_t minimumValid = static_cast<std::size_t>(
            LensingConfig::noise_plane_min_valid_fraction
            * static_cast<double>(geometricCandidates));
        if (points.size() <= std::max<std::size_t>(3U, minimumValid)) {
            return false;
        }

        UniversalUtils::findSlope2D(points, aa, bb, cc);
        return std::isfinite(aa) && std::isfinite(bb) && std::isfinite(cc);
    }

}

#endif // NOISE_PLANE_FIT_HPP
