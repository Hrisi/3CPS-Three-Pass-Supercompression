#include<iostream>
#include<fstream>
#include<algorithm>
#include<vector>
#include<sstream>

namespace nv {

    class Color32;
    struct ColorBlock;
    struct BlockDXT1;
    class Vector3;
    class Vector4;

    // All these functions return MSE.

    float compress_dxt1_single_color_optimal(Color32 c, BlockDXT1 * output);
    float compress_dxt1_single_color_optimal(const Vector3 & color, BlockDXT1 * output);

    float compress_dxt1_single_color(const Vector3 * colors, const float * weights, int count, const Vector3 & color_weights, BlockDXT1 * output);
    //float compress_dxt1_least_squares_fit(const Vector4 input_colors[16], const Vector3 * colors, const float * weights, int count, const Vector3 & color_weights, BlockDXT1 * output);
    float compress_dxt1_bounding_box_exhaustive(const Vector4 input_colors[16], const Vector3 * colors, const float * weights, int count, const Vector3 & color_weights, bool three_color_mode, int search_limit, BlockDXT1 * output);
    // --- our solution
    int division_transform(Vector3 * colors, float * weights, int count, bool isLch=false); // --not in the paper, less efficient
    void look_for_end_points_in_set(const Vector3 * unique_colors, Vector3 * colors, Vector3 * start, Vector3 * end, const float * weights, int init_count, int count, bool three_color_mode);
    // --- our solution --- end
    void compress_dxt1_cluster_fit(const Vector4 input_colors[16], const Vector3 * colors, const float * weights, int count, const Vector3 & color_weights, bool three_color_mode, BlockDXT1 * output);

    // Cluster fit end point selection.
    float compress_dxt1(const Vector4 input_colors[16], const float input_weights[16], const Vector3 & color_weights, bool three_color_mode, BlockDXT1 * output);

    // Quick end point selection followed by least squares refinement.
    float compress_dxt1_fast(const Vector4 input_colors[16], const float input_weights[16], const Vector3 & color_weights, BlockDXT1 * output);

}

