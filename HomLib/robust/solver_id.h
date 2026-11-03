#pragma once

namespace HomLib {

enum class SolverId {
    // Point-based (PoseEstimator<T>)
    FitzgibbonCVPR2001OneSided,
    FitzgibbonCVPR2001TwoSidedEqual,
    KukelovaCVPR2015TwoSidedEqual,
    KukelovaCVPR2015TwoSidedEqual6Pt,
    KukelovaCVPR2015TwoSided,
    KukelovaCVPR2015TwoSided6Pt,
    NakanoICPR2025OneSided,
    NakanoICPR2025OneSidedRight,
    Wadenback3DV2026OneSided,
    Wadenback3DV2026TwoSidedEqual,
    Wadenback3DV2026TwoSided,
    // Affine-based (AffinePoseEstimator<T>)
    BarathVISAPP2016Affine,
    ValtonenOrnhagICPR2026AffineSingleSided,
    // Orientation-based (OrientationPoseEstimator<T>)
    ValtonenOrnhagICPR2026OrientationSingleSided,
};

enum class FeatureType { Points, Affine, Orientation };

inline FeatureType feature_type_for(SolverId id) {
    switch (id) {
        case SolverId::BarathVISAPP2016Affine:
        case SolverId::ValtonenOrnhagICPR2026AffineSingleSided:
            return FeatureType::Affine;
        case SolverId::ValtonenOrnhagICPR2026OrientationSingleSided:
            return FeatureType::Orientation;
        default:
            return FeatureType::Points;
    }
}

} // namespace HomLib
