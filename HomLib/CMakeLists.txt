set(SRCS
    helpers/generate_problem_instance.cpp
    helpers/gj.cpp
    helpers/normalize2dpts.cpp
    helpers/radial.cpp
    helpers/roots.cpp
    robust/refinement.cpp
    solvers/fitzgibbon_cvpr_2001/get_fitzgibbon_cvpr_2001.cpp
    solvers/fitzgibbon_cvpr_2001/get_fitzgibbon_cvpr_2001_single.cpp
    solvers/kukelova_cvpr_2015/get_kukelova_cvpr_2015.cpp
    solvers/kukelova_cvpr_2015/get_kukelova_cvpr_2015_6pt.cpp
    solvers/kukelova_cvpr_2015/solver_kukelova_cvpr_2015.cpp
    solvers/nakano_icpr_2025/get_nakano_icpr_2025.cpp
    solvers/valtonenornhag_icpr_2020/get_valtonenornhag_icpr_2020_fHf.cpp
    solvers/valtonenornhag_icpr_2020/solver_valtonenornhag_icpr_2020_fHf.cpp
    solvers/valtonenornhag_wacv_2021/get_valtonenornhag_wacv_2021_fHf.cpp
    solvers/valtonenornhag_wacv_2021/get_valtonenornhag_wacv_2021_frHfr.cpp
    solvers/valtonenornhag_wacv_2021/solver_valtonenornhag_wacv_2021_fHf.cpp
    solvers/valtonenornhag_wacv_2021/solver_valtonenornhag_wacv_2021_frHfr.cpp
    solvers/wadenback_2025/get_wadenback_2025_double_sided.cpp
    solvers/wadenback_2025/get_wadenback_2025_double_sided_equal.cpp
    solvers/wadenback_2025/get_wadenback_2025_one_sided.cpp
    solvers/wadenback_2025/wadenback_common.cpp
)

add_library(HomLib STATIC ${SRCS})
target_link_libraries(HomLib PRIVATE Eigen3::Eigen)
target_link_libraries(HomLib PRIVATE PoseLib::PoseLib)
target_include_directories(HomLib PRIVATE ${HOMLIB_INCLUDE_DIRS})
target_compile_options(HomLib PRIVATE -fPIC)