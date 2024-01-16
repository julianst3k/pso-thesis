#pragma once
#include <vector>
#include <cmath>
#include <algorithm>
#include "parameter_aggregate.h"
#include <pybind11/numpy.h>
#include <random>
#include <memory>

namespace shadowing{
struct Width_Height_Prob{
    float width;
    float height;
    float probability;
    float ymin;
    float ymax;
    float xmax;
    float ymax_real;
    float xmin;
};


struct WH_Probabilities{
    std::vector<Width_Height_Prob*> WH_collection;
    void push_probability(float width, float height, float probability, float ymin, float ymax, float ymax_real, float xmax, float xmin){
        Width_Height_Prob* whp = new Width_Height_Prob();
        whp->width = width; whp->height = height; whp->probability = probability;
        whp->ymax_real = ymax_real; whp->xmin = xmin;
        whp->ymin = ymin; whp->ymax = ymax; whp->xmax = xmax;
        this->WH_collection.push_back(whp);
    };
};


struct Shadowing_Parameters{
    float width;
    float height;
    float ymin;
    float xmax;
    float xmin;
    float ymax;
    float zj;
    float zi;
    float probability;
    float ymax_real;
    float ymin_real;
    std::vector<float> coordinate_transmitter;
    std::vector<float> coordinate_receiver;
    Shadowing_Parameters(std::vector<float> transmitter, std::vector<float> receiv, float width, float height, float ymin, float ymin_real, float ymax, float ymax_real, float xmax, float xmin, float zi, float zj, float probability):
    coordinate_transmitter(transmitter), coordinate_receiver(receiv), width(width), height(height), ymin(ymin), ymax(ymax), ymin_real(ymin_real), ymax_real(ymax_real), xmax(xmax), xmin(xmin), zi(zi), zj(zj), probability(probability) {};
    Shadowing_Parameters(float xt, float yt, float xr, float yr, float width, float height, float ymin, float ymax, float xmax, float xmin, float ymin_real, float ymax_real, float zi, float zj, float probability):
    coordinate_transmitter(std::vector{xt,yt,zi}), coordinate_receiver(std::vector{xr,yr,zj}), width(width), height(height), ymin(ymin), ymin_real(ymin_real), ymax(ymax), ymax_real(ymax_real), xmax(xmax), xmin(xmin), zi(zi), zj(zj), probability(probability) {};
    void print_stuff(){

    };
    
};

struct Shadowing_Secondary_Parameters{
    float W1;
    float W2;
    float H1;
    float H2;
    float alpha;
    float beta;
    float ymax;
    float ymin;
    float ymax_real;
    float sind;
    float cotsind;

    float xmax;
    void calculate_parameters(Shadowing_Parameters* parameters){
        float x_i = parameters->coordinate_transmitter[0];
        float x_j = parameters->coordinate_receiver[0];
        float y_i = parameters->coordinate_transmitter[1];
        float y_j = parameters->coordinate_receiver[1];
        float z_i = parameters->coordinate_transmitter[2];
        float z_j = parameters->coordinate_receiver[2];
        float x_ij = x_i - x_j;
        if(x_ij == 0){
            x_ij = x_ij + 0.00001;
            x_i = x_i + 0.00001;
        };
        float y_ij = y_i - y_j;
        if(y_ij == 0){
            y_ij = y_ij - 0.00001;
            y_i = y_i - 0.00001;
        };
        float z_ij = z_i - z_j;
        float d_ij = std::sqrt(std::pow((x_i-x_j),2)+std::pow((y_i-y_j),2)+std::pow(z_i-z_j,2));
        float d_ij_nz = std::sqrt(std::pow((x_i-x_j),2)+std::pow((y_i-y_j),2));
        this->sind = calculate_sind(x_ij, y_ij, z_ij);
        this->cotsind = calculate_cotsind(x_ij, y_ij, z_ij);
        this->W1 = calculate_W1(x_i, x_j, y_i, y_j, x_ij);
        this->W2 = calculate_W2(d_ij_nz, parameters->width, x_ij);
        this->alpha = calculate_alpha(x_ij, y_ij);
        this->beta = 1 + 1/std::pow(this->alpha,2);
        this->H1 = calculate_H1(parameters->height, parameters->zj, d_ij, y_ij, this->beta)*this->cotsind;
        this->H2 = calculate_H2(y_ij, x_ij, x_j, y_j, x_i, y_i, this->sind, d_ij, this->beta);
        this->ymax = parameters->ymax;
        this->ymin = parameters->ymin;
        this->xmax = parameters->xmax;
        this->ymax_real = parameters->ymax_real;
    };
    void print_stuff(){
        std::cout << "W1/alpha " << std::to_string(this->W1/this->alpha) << "\n";
        std::cout << "W2/alpha " << std::to_string(this->W2/this->alpha) << "\n"; 
        std::cout << "W2/alpha + W1/alpha" << std::to_string(this->W2/this->alpha+this->W1/this->alpha) << "\n"; 
        std::cout << "W2/alpha - W1/alpha" << std::to_string(this->W2/this->alpha-this->W1/this->alpha) << "\n"; 
        // std::cout << "H1 " << std::to_string(this->H1) << "\n";
     // std::cout << "H2 " << std::to_string(this->H2) << "\n";
       // std::cout << "Alph " << std::to_string(this->alpha) << "\n";
        // std::cout << "Beta " << std::to_string(this->beta) << "\n";
        // std::cout << "cotsind " << std::to_string(this->cotsind) << "\n";
        //std::cout << "sind " << std::to_string(this->sind) << "\n";
    };
    float calculate_sind(float x_ij, float y_ij, float z_ij);
    float calculate_cotsind(float x_ij, float y_ij, float z_ij);
    float calculate_W1(float x_i, float x_j, float y_i, float y_j, float x_ij);
    float calculate_W2(float d_ij, float w, float x_ij);
    float calculate_H1(float h, float zj, float d_ij, float y_ij, float beta);
    float calculate_H2(float y_ij, float x_ij, float x_j, float y_j, float x_i, float y_i, float sind, float d_ij, float beta);
    float calculate_alpha(float x_ij, float y_ij);
    
};
struct Shadowing_Parameters_Coll{
    std::vector<std::unique_ptr<Shadowing_Parameters>> SP_collection;
    std::vector<std::unique_ptr<Shadowing_Secondary_Parameters>> SPS_collection;
    void create_collection(float xt, float yt, float zi, float xr, float yr, float zj, WH_Probabilities* wh_coll){
        float min_z; float max_z; float ymax;
        if(zi<zj){
            min_z = zi;
            max_z = zj;
        }else{
            min_z = zj;
            max_z = zi;
        }
        if(yt<yr){
            ymax = yr;
        }
        else{
            ymax = yt;
        }

        for(Width_Height_Prob* whp : wh_coll->WH_collection){
            if(yr<yt){
                std::unique_ptr<Shadowing_Parameters> shadow =  std::make_unique<Shadowing_Parameters>(xr, yr, xt, yt, whp->width, whp->height, yr,  ymax, whp->xmax, whp->xmin, whp->ymin, whp->ymax_real, max_z, min_z, whp->probability);
                this->SPS_collection.emplace_back(std::move(shadowing_param(shadow.get())));
                this->SP_collection.emplace_back(std::move(shadow));

            }else{
                std::unique_ptr<Shadowing_Parameters> shadow =  std::make_unique<Shadowing_Parameters>(xt, yt, xr, yr, whp->width, whp->height, yt, ymax, whp->xmax, whp->xmin, whp->ymin, whp->ymax_real, max_z, min_z, whp->probability);
                this->SPS_collection.emplace_back(std::move(shadowing_param(shadow.get())));
                this->SP_collection.emplace_back(std::move(shadow));

            };
        };
    };
    std::unique_ptr<Shadowing_Secondary_Parameters> shadowing_param(Shadowing_Parameters* shadow){
        std::unique_ptr<Shadowing_Secondary_Parameters> shadow_param = std::make_unique<Shadowing_Secondary_Parameters>();
        shadow_param->calculate_parameters(shadow);
        return shadow_param;
    };
};
float calculate_expectancy(shadowing::Shadowing_Parameters_Coll* shadow_param);

};
float calculate_total_integral(shadowing::Shadowing_Parameters* parameters, shadowing::Shadowing_Secondary_Parameters* shadow_secondary);
float calculate_total_integral(shadowing::Shadowing_Parameters* parameters, shadowing::Shadowing_Secondary_Parameters* shadow_secondary);
bool shadowing_comparator(std::vector<float>* value_set, float anchor);
float max_set(std::vector<float>* value_set);
float min_set(std::vector<float>* value_set);


enum Sign_Order {LEQ, GEQ};