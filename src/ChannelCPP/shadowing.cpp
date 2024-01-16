#include <vector>
#include <cmath>
#include <algorithm>
#include <pybind11/numpy.h>
#include <random>
#include <cfloat>
#include "shadowing.h"
#include <memory>

namespace sh = shadowing;



float sh::Shadowing_Secondary_Parameters::calculate_W1(float x_i, float x_j, float y_i, float y_j, float x_ij){

    return (x_i*y_j-x_j*y_i)/x_ij;
}

float sh::Shadowing_Secondary_Parameters::calculate_W2(float d_ij, float w, float x_ij){
    return (w*d_ij)/(2*x_ij);
}

float sh::Shadowing_Secondary_Parameters::calculate_H1(float h, float zj, float d_ij, float y_ij, float beta){
    return (h-zj)*d_ij/(y_ij*beta);
}

float sh::Shadowing_Secondary_Parameters::calculate_H2(float y_ij, float x_ij, float x_j, float y_j, float x_i, float y_i, float sind, float d_ij, float beta){
    return (std::pow(d_ij*sind,2)+std::pow(y_j,2)*beta-std::pow(y_i,2)*beta)/(2*y_ij*beta);
}

float sh::Shadowing_Secondary_Parameters::calculate_alpha(float x_ij, float y_ij){
    return y_ij/x_ij;
}
float sh::Shadowing_Secondary_Parameters::calculate_sind(float x_ij, float y_ij, float z_ij){
    float d = std::sqrt(std::pow(x_ij,2)+std::pow(y_ij,2)+std::pow(z_ij,2));
    float x = z_ij/d;
    return std::sqrt(1-std::pow(x,2));
}
float sh::Shadowing_Secondary_Parameters::calculate_cotsind(float x_ij, float y_ij, float z_ij){
    float d = std::sqrt(std::pow(x_ij,2)+std::pow(y_ij,2)+std::pow(z_ij,2));
    float x = z_ij/d;
    return (1-std::pow(x,2))/x;
}



bool shadowing_comparator(std::vector<float>* value_set, float anchor, Sign_Order sign){
    bool ret = true;
    for(float value: *value_set){
        if(sign == GEQ){
            ret = ret && value >= anchor;
        } else{
            ret = ret && value <= anchor;
        }
    }
    return ret;
}

float max_set(std::vector<float>* value_set){
    float max = -FLT_MAX;
    for(float value: *value_set){
        if(value > max){
            max = value;
        } 
    }
    return max;

}
float min_set(std::vector<float>* value_set){
    float min = FLT_MAX;
    for(float value: *value_set){
        if(value < min){
            min = value;
        } 
    }
    return min;

}
float quadratic_function(float a, float b, float x){
    return a/2*std::pow(x,2)+b*x;  
}
float quadratic_function_diff(float a, float b, float lb, float ub){
    return quadratic_function(a,b,ub)-quadratic_function(a,b,lb);
}
float sh::calculate_expectancy(sh::Shadowing_Parameters_Coll* shadow_param){
    float tot = 0;
    int i = 0;
    for(std::unique_ptr<sh::Shadowing_Parameters>& parameters : shadow_param -> SP_collection){
        float norm = (parameters->xmax-parameters->xmin)*(parameters->ymax_real-parameters->ymin_real);
        tot += calculate_total_integral(parameters.get(), (shadow_param->SPS_collection[i]).get())*parameters->probability/norm;
        i += 1;
    }
    return tot;
}

void print_vector(std::vector<float>* vector, const std::string& val){
    std::cout << val << "\n";
    for(float e: *vector){
        std::cout << std::to_string(e) << "\n";
     }
}

float integral_11(sh::Shadowing_Parameters* parameters, sh::Shadowing_Secondary_Parameters* shadow_parameters){
    float ymax = parameters->ymax;
    float xmax = parameters->xmax;
    float minv = parameters->ymin;
    float xmin = parameters->xmin;

    if(minv<shadow_parameters->H1-shadow_parameters->H2){
        minv = shadow_parameters->H1-shadow_parameters->H2;
    }
    std::vector<float> min_vector{parameters->xmax, (minv-shadow_parameters->W1)/shadow_parameters->alpha, 
    (parameters->ymax-shadow_parameters->W1)/shadow_parameters->alpha};
    std::vector<float> max_vector{xmin, (shadow_parameters->ymax+shadow_parameters->W2-shadow_parameters->W1)/shadow_parameters->alpha};
    float diff = quadratic_function_diff(0, ymax, max_set(&max_vector), min_set(&min_vector));

    float total = diff*shadowing_comparator(&min_vector, xmin, GEQ)*shadowing_comparator(&max_vector, xmax, LEQ);

    return total;
}

float integral_12(sh::Shadowing_Parameters* parameters, sh::Shadowing_Secondary_Parameters* shadow_parameters){
    float ymax = parameters->ymax;
    float xmax = parameters->xmax;
    float minv = parameters->ymin;
    float xmin = parameters->xmin;

    if(minv<shadow_parameters->H1-shadow_parameters->H2){
        minv = shadow_parameters->H1-shadow_parameters->H2;
    }
    std::vector<float> min_vector{parameters->xmax, (minv-shadow_parameters->W1)/shadow_parameters->alpha};
    std::vector<float> max_vector{xmin, (shadow_parameters->ymax+shadow_parameters->W2-shadow_parameters->W1)/shadow_parameters->alpha, 
    (parameters->ymax-shadow_parameters->W1)/shadow_parameters->alpha};
    float diff = quadratic_function_diff(shadow_parameters->alpha, shadow_parameters->W1, max_set(&max_vector), min_set(&min_vector));

    float total = diff*shadowing_comparator(&min_vector, xmin, GEQ)*shadowing_comparator(&max_vector, xmax, LEQ);

    return total;
}

float integral_13(sh::Shadowing_Parameters* parameters, sh::Shadowing_Secondary_Parameters* shadow_parameters){
    float ymin = parameters->ymin;
    float xmax = parameters->xmax;
    float minv = parameters->ymin;
    float xmin = parameters->xmin;

    if(minv<shadow_parameters->H1-shadow_parameters->H2){
        minv = shadow_parameters->H1-shadow_parameters->H2;
    }
    std::vector<float> min_vector{parameters->xmax, (minv-shadow_parameters->W1)/shadow_parameters->alpha};
    std::vector<float> max_vector{xmin, (shadow_parameters->ymax+shadow_parameters->W2-shadow_parameters->W1)/shadow_parameters->alpha, 
    (minv+shadow_parameters->W2-shadow_parameters->W1)/shadow_parameters->alpha};
    float diff = quadratic_function_diff(0, minv, max_set(&max_vector), min_set(&min_vector));
    float total = diff*shadowing_comparator(&min_vector, xmin, GEQ)*shadowing_comparator(&max_vector, xmax, LEQ);

    return total;
}

float integral_14(sh::Shadowing_Parameters* parameters, sh::Shadowing_Secondary_Parameters* shadow_parameters){
    float ymax = parameters->ymax;
    float xmax = parameters->xmax;
    float minv = parameters->ymin;
    float xmin = parameters->xmin;

    if(minv<shadow_parameters->H1-shadow_parameters->H2){
        minv = shadow_parameters->H1-shadow_parameters->H2;
    }
    std::vector<float> min_vector{parameters->xmax, (minv-shadow_parameters->W1)/shadow_parameters->alpha,
    (minv+shadow_parameters->W2-shadow_parameters->W1)/shadow_parameters->alpha};
    std::vector<float> max_vector{xmin, (shadow_parameters->ymax+shadow_parameters->W2-shadow_parameters->W1)/shadow_parameters->alpha};
    float diff = quadratic_function_diff(shadow_parameters->alpha, shadow_parameters->W1-shadow_parameters->W2, max_set(&max_vector), min_set(&min_vector));
    float total = diff*shadowing_comparator(&min_vector, xmin, GEQ)*shadowing_comparator(&max_vector, xmax, LEQ);

    return total;
}

float integral_21(sh::Shadowing_Parameters* parameters, sh::Shadowing_Secondary_Parameters* shadow_parameters){
    float ymax = parameters->ymax;
    float xmax = parameters->xmax;
    float minv = parameters->ymin;
    float xmin = parameters->xmin;

    if(minv<shadow_parameters->H1-shadow_parameters->H2){
        minv = shadow_parameters->H1-shadow_parameters->H2;
    }
    std::vector<float> min_vector{parameters->xmax, (parameters->ymax-shadow_parameters->W1)/shadow_parameters->alpha};
    std::vector<float> max_vector{xmin, (minv+shadow_parameters->W2-shadow_parameters->W1)/shadow_parameters->alpha,
    (shadow_parameters->ymax+shadow_parameters->W2-shadow_parameters->W1)/shadow_parameters->alpha};
    float diff = quadratic_function_diff(0, ymax, max_set(&max_vector), min_set(&min_vector));
    float total = diff*shadowing_comparator(&min_vector, xmin, GEQ)*shadowing_comparator(&max_vector, xmax, LEQ);

    return total;
}


float integral_22(sh::Shadowing_Parameters* parameters, sh::Shadowing_Secondary_Parameters* shadow_parameters){
    float ymax = parameters->ymax;
    float xmax = parameters->xmax;
    float minv = parameters->ymin;
    float xmin = parameters->xmin;

    if(minv<shadow_parameters->H1-shadow_parameters->H2){
        minv = shadow_parameters->H1-shadow_parameters->H2;
    }    std::vector<float> min_vector{parameters->xmax, (parameters->ymax-shadow_parameters->W1)/shadow_parameters->alpha,
    (shadow_parameters->ymax+shadow_parameters->W2-shadow_parameters->W1)/shadow_parameters->alpha};
    std::vector<float> max_vector{xmin, (minv+shadow_parameters->W2-shadow_parameters->W1)/shadow_parameters->alpha};
    float diff = quadratic_function_diff(shadow_parameters->alpha, shadow_parameters->W1-shadow_parameters->W2, max_set(&max_vector), min_set(&min_vector));
    float total = diff*shadowing_comparator(&min_vector, xmin, GEQ)*shadowing_comparator(&max_vector, xmax, LEQ);

    return total;
}

float integral_23(sh::Shadowing_Parameters* parameters, sh::Shadowing_Secondary_Parameters* shadow_parameters){
    float ymin = parameters->ymin;
    float xmax = parameters->xmax;
    float minv = parameters->ymin;
    float xmin = parameters->xmin;
    if(minv<shadow_parameters->H1-shadow_parameters->H2){
        minv = shadow_parameters->H1-shadow_parameters->H2;
    }
    std::vector<float> min_vector{parameters->xmax, (parameters->ymax-shadow_parameters->W1)/shadow_parameters->alpha,
    (minv-shadow_parameters->W1)/shadow_parameters->alpha};
    std::vector<float> max_vector{xmin, (minv+shadow_parameters->W2-shadow_parameters->W1)/shadow_parameters->alpha};
    float diff = quadratic_function_diff(0, minv, max_set(&max_vector), min_set(&min_vector));
    float total = diff*shadowing_comparator(&min_vector, xmin, GEQ)*shadowing_comparator(&max_vector, xmax, LEQ);

    return total;
}

float integral_24(sh::Shadowing_Parameters* parameters, sh::Shadowing_Secondary_Parameters* shadow_parameters){
    float ymin = parameters->ymin;
    float xmax = parameters->xmax;
    float xmin = parameters->xmin;
    float minv = parameters->ymin;
    if(minv<shadow_parameters->H1-shadow_parameters->H2){
        minv = shadow_parameters->H1-shadow_parameters->H2;
    }
    std::vector<float> min_vector{parameters->xmax, (parameters->ymax-shadow_parameters->W1)/shadow_parameters->alpha};
    std::vector<float> max_vector{parameters->xmin, (minv+shadow_parameters->W2-shadow_parameters->W1)/shadow_parameters->alpha,
    (minv-shadow_parameters->W1)/shadow_parameters->alpha};
    float diff = quadratic_function_diff(shadow_parameters->alpha, shadow_parameters->W1, max_set(&max_vector), min_set(&min_vector));
    float total = diff*shadowing_comparator(&min_vector, xmin, GEQ)*shadowing_comparator(&max_vector, xmax, LEQ);

    return total;
}


float integral_31(sh::Shadowing_Parameters* parameters, sh::Shadowing_Secondary_Parameters* shadow_parameters){
    float ymax = parameters->ymax;
    float xmax = parameters->xmax;
    float minv = parameters->ymin;
    float xmin = parameters->xmin;

    if(minv<shadow_parameters->H1-shadow_parameters->H2){
        minv = shadow_parameters->H1-shadow_parameters->H2;
    }
    std::vector<float> min_vector{parameters->xmax, (parameters->ymax-shadow_parameters->W2-shadow_parameters->W1)/shadow_parameters->alpha,
    (minv-shadow_parameters->W2-shadow_parameters->W1)/shadow_parameters->alpha};
    std::vector<float> max_vector{xmin, (shadow_parameters->ymax-shadow_parameters->W1)/shadow_parameters->alpha};
    float diff = quadratic_function_diff(0, ymax, max_set(&max_vector), min_set(&min_vector));


    float total = diff*shadowing_comparator(&min_vector, xmin, GEQ)*shadowing_comparator(&max_vector, xmax, LEQ);

    return total;
}

float integral_32(sh::Shadowing_Parameters* parameters, sh::Shadowing_Secondary_Parameters* shadow_parameters){
    float ymax = parameters->ymax;
    float xmax = parameters->xmax;
    float minv = parameters->ymin;
    float xmin = parameters->xmin;

    if(minv<shadow_parameters->H1-shadow_parameters->H2){
        minv = shadow_parameters->H1-shadow_parameters->H2;
    }
    std::vector<float> min_vector{parameters->xmax,
    (minv-shadow_parameters->W2-shadow_parameters->W1)/shadow_parameters->alpha};
    std::vector<float> max_vector{xmin, (shadow_parameters->ymax-shadow_parameters->W1)/shadow_parameters->alpha,
    (parameters->ymax-shadow_parameters->W2-shadow_parameters->W1)/shadow_parameters->alpha};
    float diff = quadratic_function_diff(shadow_parameters->alpha, shadow_parameters->W1+shadow_parameters->W2, max_set(&max_vector), min_set(&min_vector));


    float total = diff*shadowing_comparator(&min_vector, xmin, GEQ)*shadowing_comparator(&max_vector, xmax, LEQ);

    return total;
}


float integral_33(sh::Shadowing_Parameters* parameters, sh::Shadowing_Secondary_Parameters* shadow_parameters){
    float ymin = parameters->ymin;
    float xmax = parameters->xmax;
    float minv = parameters->ymin;
    float xmin = parameters->xmin;
    if(minv<shadow_parameters->H1-shadow_parameters->H2){
        minv = shadow_parameters->H1-shadow_parameters->H2;
    }
    std::vector<float> min_vector{parameters->xmax,
    (minv-shadow_parameters->W2-shadow_parameters->W1)/shadow_parameters->alpha};
    std::vector<float> max_vector{xmin, (shadow_parameters->ymax-shadow_parameters->W1)/shadow_parameters->alpha,
    (minv-shadow_parameters->W1)/shadow_parameters->alpha};
    float diff = quadratic_function_diff(0, minv, max_set(&max_vector), min_set(&min_vector));
    float total = diff*shadowing_comparator(&min_vector, xmin, GEQ)*shadowing_comparator(&max_vector, xmax, LEQ);

    return total;
}


float integral_34(sh::Shadowing_Parameters* parameters, sh::Shadowing_Secondary_Parameters* shadow_parameters){
    float ymin = parameters->ymin;
    float xmax = parameters->xmax;
    float minv = parameters->ymin;
    float xmin = parameters->xmin;
    if(minv<shadow_parameters->H1-shadow_parameters->H2){
        minv = shadow_parameters->H1-shadow_parameters->H2;
    }
    std::vector<float> min_vector{parameters->xmax,
    (minv-shadow_parameters->W2-shadow_parameters->W1)/shadow_parameters->alpha,
    (minv-shadow_parameters->W1)/shadow_parameters->alpha};
    std::vector<float> max_vector{xmin, (shadow_parameters->ymax-shadow_parameters->W1)/shadow_parameters->alpha};
    float diff = quadratic_function_diff(shadow_parameters->alpha, shadow_parameters->W1, max_set(&max_vector), min_set(&min_vector));
    float total = diff*shadowing_comparator(&min_vector, xmin, GEQ)*shadowing_comparator(&max_vector, xmax, LEQ);
    return total;
}


float integral_41(sh::Shadowing_Parameters* parameters, sh::Shadowing_Secondary_Parameters* shadow_parameters){
    float ymax = parameters->ymax;
    float xmax = parameters->xmax;
    float minv = parameters->ymin;
    float xmin = parameters->xmin;

    if(minv<shadow_parameters->H1-shadow_parameters->H2){
        minv = shadow_parameters->H1-shadow_parameters->H2;
    }
    std::vector<float> min_vector{parameters->xmax, (parameters->ymax-shadow_parameters->W2-shadow_parameters->W1)/shadow_parameters->alpha};
    std::vector<float> max_vector{xmin, (shadow_parameters->ymax-shadow_parameters->W1)/shadow_parameters->alpha,
    (minv-shadow_parameters->W1)/shadow_parameters->alpha};
    float diff = quadratic_function_diff(0, ymax, max_set(&max_vector), min_set(&min_vector));
    float total = diff*shadowing_comparator(&min_vector, xmin, GEQ)*shadowing_comparator(&max_vector, xmax, LEQ);
    
    return total;
}

float integral_42(sh::Shadowing_Parameters* parameters, sh::Shadowing_Secondary_Parameters* shadow_parameters){
    float ymax = parameters->ymax;
    float xmax = parameters->xmax;
    float minv = parameters->ymin;
    float xmin = parameters->xmin;

    if(minv<shadow_parameters->H1-shadow_parameters->H2){
        minv = shadow_parameters->H1-shadow_parameters->H2;
    }    
    std::vector<float> min_vector{parameters->xmax, (parameters->ymax-shadow_parameters->W2-shadow_parameters->W1)/shadow_parameters->alpha,
    (shadow_parameters->ymax-shadow_parameters->W1)/shadow_parameters->alpha};
    std::vector<float> max_vector{xmin, (minv-shadow_parameters->W1)/shadow_parameters->alpha};
    float diff = quadratic_function_diff(shadow_parameters->alpha, shadow_parameters->W1, max_set(&max_vector), min_set(&min_vector));
    float total = diff*shadowing_comparator(&min_vector, xmin, GEQ)*shadowing_comparator(&max_vector, xmax, LEQ);

    return total;
}

float integral_43(sh::Shadowing_Parameters* parameters, sh::Shadowing_Secondary_Parameters* shadow_parameters){
    float ymin = parameters->ymin;
    float xmax = parameters->xmax;
    float minv = parameters->ymin;
    float xmin = parameters->xmin;

    if(minv<shadow_parameters->H1-shadow_parameters->H2){
        minv = shadow_parameters->H1-shadow_parameters->H2;
    }
    std::vector<float> min_vector{parameters->xmax,
    (minv-shadow_parameters->W2-shadow_parameters->W1)/shadow_parameters->alpha,
    (parameters->ymax-shadow_parameters->W2-shadow_parameters->W1)/shadow_parameters->alpha};
    std::vector<float> max_vector{xmin, (minv-shadow_parameters->W1)/shadow_parameters->alpha};
    float diff = quadratic_function_diff(0, minv, max_set(&max_vector), min_set(&min_vector));
    float total = diff*shadowing_comparator(&min_vector, xmin, GEQ)*shadowing_comparator(&max_vector, xmax, LEQ);

    return total;
}

float integral_44(sh::Shadowing_Parameters* parameters, sh::Shadowing_Secondary_Parameters* shadow_parameters){
    float ymin = parameters->ymin;
    float xmax = parameters->xmax;
    float minv = parameters->ymin;
    float xmin = parameters->xmin;

    if(minv<shadow_parameters->H1-shadow_parameters->H2){
        minv = shadow_parameters->H1-shadow_parameters->H2;
    }
    std::vector<float> min_vector{parameters->xmax,
    (parameters->ymax-shadow_parameters->W2-shadow_parameters->W1)/shadow_parameters->alpha};
    std::vector<float> max_vector{xmin, (minv-shadow_parameters->W1)/shadow_parameters->alpha,
    (minv-shadow_parameters->W2-shadow_parameters->W1)/shadow_parameters->alpha};
    float diff = quadratic_function_diff(shadow_parameters->alpha, shadow_parameters->W2+shadow_parameters->W1, max_set(&max_vector), min_set(&min_vector));
    float total = diff*shadowing_comparator(&min_vector, xmin, GEQ)*shadowing_comparator(&max_vector, xmax, LEQ);

    return total;
}



float calculate_total_integral_13(sh::Shadowing_Parameters* parameters, sh::Shadowing_Secondary_Parameters* shadow_secondary){

    float integral_one = integral_11(parameters, shadow_secondary)+integral_12(parameters, shadow_secondary)-(integral_13(parameters, shadow_secondary)+
    integral_14(parameters, shadow_secondary));
    float integral_three = integral_31(parameters, shadow_secondary)+integral_32(parameters, shadow_secondary)-(integral_33(parameters, shadow_secondary)+
    integral_34(parameters, shadow_secondary));
    float ret = integral_one + integral_three;
    return ret;
}
float calculate_total_integral_24(sh::Shadowing_Parameters* parameters, sh::Shadowing_Secondary_Parameters* shadow_secondary){

    float integral_two = integral_21(parameters, shadow_secondary)+integral_22(parameters, shadow_secondary)-(integral_23(parameters, shadow_secondary)+
    integral_24(parameters, shadow_secondary));
    float integral_four = integral_41(parameters, shadow_secondary)+integral_42(parameters, shadow_secondary)-(integral_43(parameters, shadow_secondary)+
    integral_44(parameters, shadow_secondary));
    
    return integral_two + integral_four;
}

float calculate_total_integral(sh::Shadowing_Parameters* parameters, sh::Shadowing_Secondary_Parameters* shadow_secondary){
        float x_i = parameters->coordinate_transmitter[0];
        float x_j = parameters->coordinate_receiver[0]; 
        float x_ij = x_i - x_j;

        if(shadow_secondary->H1-shadow_secondary->H2<=shadow_secondary->ymax){
            if(x_ij >= 0){
                return std::abs(calculate_total_integral_13(parameters, shadow_secondary));
            } else {
                return std::abs(calculate_total_integral_24(parameters, shadow_secondary));
            }
        }
        else{
            return 0;
        }
}
