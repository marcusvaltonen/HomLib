#pragma once

#include <stdexcept>

#include <solver_id.h>
#include <get_fitzgibbon_cvpr_2001.hpp>
#include <get_kukelova_cvpr_2015.hpp>
#include <get_nakano_icpr_2025.hpp>
#include <get_wadenback_3dv_2026.hpp>
#include <get_barath_visapp_2016.hpp>
#include <get_valtonenornhag_icpr_2026.hpp>

namespace HomLib {

template <typename Func>
auto dispatch_point_solver(SolverId id, Func&& func) {
    switch (id) {
        case SolverId::FitzgibbonCVPR2001OneSided:
            return func(FitzgibbonCVPR2001::SolverSingleSided{});
        case SolverId::FitzgibbonCVPR2001TwoSidedEqual:
            return func(FitzgibbonCVPR2001::SolverTwoSidedEqual{});
        case SolverId::KukelovaCVPR2015TwoSidedEqual:
            return func(KukelovaCVPR2015::SolverTwoSidedEqual{});
        case SolverId::KukelovaCVPR2015TwoSidedEqual6Pt:
            return func(KukelovaCVPR2015::SolverTwoSidedEqual6Pt{});
        case SolverId::KukelovaCVPR2015TwoSided:
            return func(KukelovaCVPR2015::SolverTwoSided{});
        case SolverId::KukelovaCVPR2015TwoSided6Pt:
            return func(KukelovaCVPR2015::SolverTwoSided6Pt{});
        case SolverId::NakanoICPR2025OneSided:
            return func(NakanoICPR2025::SolverSingleSided{});
        case SolverId::NakanoICPR2025OneSidedRight:
            return func(NakanoICPR2025::SolverSingleSidedRight{});
        case SolverId::Wadenback3DV2026OneSided:
            return func(Wadenback3DV2026::SolverSingleSided{});
        case SolverId::Wadenback3DV2026TwoSidedEqual:
            return func(Wadenback3DV2026::SolverTwoSidedEqual{});
        case SolverId::Wadenback3DV2026TwoSided:
            return func(Wadenback3DV2026::SolverTwoSided{});
        default:
            throw std::invalid_argument("Not a point-based solver");
    }
}

template <typename Func>
auto dispatch_affine_solver(SolverId id, Func&& func) {
    switch (id) {
        case SolverId::BarathVISAPP2016Affine:
            return func(BarathVISAPP2016::AffineSolver{});
        case SolverId::ValtonenOrnhagICPR2026AffineSingleSided:
            return func(ValtonenOrnhagICPR2026::AffineSolverSingleSided{});
        default:
            throw std::invalid_argument("Not an affine-based solver");
    }
}

template <typename Func>
auto dispatch_orientation_solver(SolverId id, Func&& func) {
    switch (id) {
        case SolverId::ValtonenOrnhagICPR2026OrientationSingleSided:
            return func(ValtonenOrnhagICPR2026::OrientationSolverSingleSided{});
        default:
            throw std::invalid_argument("Not an orientation-based solver");
    }
}

template <typename Func>
auto dispatch_solver(SolverId id, Func&& func) {
    switch (feature_type_for(id)) {
        case FeatureType::Points:
            return dispatch_point_solver(id, std::forward<Func>(func));
        case FeatureType::Affine:
            return dispatch_affine_solver(id, std::forward<Func>(func));
        case FeatureType::Orientation:
            return dispatch_orientation_solver(id, std::forward<Func>(func));
        default:
            throw std::invalid_argument("Unknown feature type");
    }
}

} // namespace HomLib
