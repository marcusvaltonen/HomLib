#pragma once

namespace HomLib {

enum DistortionCase {
    NO_DISTORTION,
    ONE_SIDED_LEFT,
    ONE_SIDED_RIGHT,
    TWO_SIDED_EQUAL,
    TWO_SIDED
};

enum class SolverId {
    // Point-based (PoseEstimator<T>)
    CaiPAMI2025,
    DLT,
    FitzgibbonCVPR2001OneSided,
    FitzgibbonCVPR2001TwoSidedEqual,
    KukelovaCVPR2015TwoSidedEqual,
    KukelovaCVPR2015TwoSidedEqual6Pt,
    KukelovaCVPR2015TwoSided,
    KukelovaCVPR2015TwoSided6Pt,
    NakanoICPR2025OneSided,
    NakanoICPR2025OneSidedRight,
    FitzgibbonCVPR2001OneSidedRight,
    Wadenback3DV2026OneSided,
    Wadenback3DV2026OneSidedRight,
    Wadenback3DV2026TwoSidedEqual,
    Wadenback3DV2026TwoSided,
    GuoPRL2019,
    // Affine-based (AffinePoseEstimator<T>)
    BarathVISAPP2016Affine,
    ValtonenOrnhagICPR2026AffineOneSided,
    ValtonenOrnhagICPR2026AffineOneSidedRight,
    // Orientation-based (OrientationPoseEstimator<T>)
    ValtonenOrnhagICPR2026OrientationOneSided,
    ValtonenOrnhagICPR2026OrientationOneSidedRight,
};

enum class FeatureType { Points, Affine, Orientation };

inline FeatureType feature_type_for(SolverId id) {
    switch (id) {
        case SolverId::BarathVISAPP2016Affine:
        case SolverId::ValtonenOrnhagICPR2026AffineOneSided:
        case SolverId::ValtonenOrnhagICPR2026AffineOneSidedRight:
            return FeatureType::Affine;
        case SolverId::ValtonenOrnhagICPR2026OrientationOneSided:
        case SolverId::ValtonenOrnhagICPR2026OrientationOneSidedRight:
            return FeatureType::Orientation;
        default:
            return FeatureType::Points;
    }
}

} // namespace HomLib
