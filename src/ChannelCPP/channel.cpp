#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <vector>
#include <cmath>
#include <algorithm>
#include <ctime>
#include <memory>
#include "parameter_aggregate.h"
#include <random>
#include "channel.h"
#include <iostream>
#include <thread>
#include "shadowing.h"
#define SIZE 8
#define pi 3.1415926535
#define RESP_SIZE 180
#define INTENSITY 5
#define TIME 0.1

namespace sh = shadowing;
namespace pa = parameter_aggregate;
namespace py = pybind11;
std::default_random_engine generator;
std::uniform_real_distribution<float> distribution(-90,90);

struct final_response_wrapper{
    std::unique_ptr<std::vector<float>> final_response;
};

void initialize_3d_matrix(float*** out_matrix, int size, int h, int w){
    out_matrix = new float**[size];
    for(int x=0; x<size; x++){
        out_matrix[x] = new float*[h];
        for(int y=0; y<h; y++){
            out_matrix[x][y] = new float[w];
        }
    }
}

enum ResponseMode {full, scattering, los, nlos};

float led_pd_channel(pa::WallParameters* wall, pa::TransmitterParameters* trans, pa::ReceiverParameters* recv, pa::TunnelParameters* tunnel,
pa::SimulationParameters* simulation, sh::WH_Probabilities* wh_coll){
    std::vector<float> final_response(simulation->time.size()+simulation->h_led.size()-1);
    HLos_Vector(wall, trans, recv, tunnel, simulation, &final_response, wh_coll);
    HNLos_Vector(wall, trans, recv, tunnel, simulation, &final_response, 0.2, wh_coll);
    HNLos_Vector(wall, trans, recv, tunnel, simulation, &final_response, 2.8, wh_coll);
    HScatter_Vector(wall, trans, recv, tunnel, simulation, &final_response, 2.8, wh_coll);
    float sum = std::reduce(final_response.begin(), final_response.end());

    final_response = std::vector<float>();
    return sum;
}


std::vector<float> led_pd_response(pa::WallParameters* wall, pa::TransmitterParameters* trans, pa::ReceiverParameters* recv, pa::TunnelParameters* tunnel,
pa::SimulationParameters* simulation, sh::WH_Probabilities* wh_coll, ResponseMode mode){
    std::vector<float> final_response(simulation->time.size()+simulation->h_led.size()-1);
    switch(mode){
        case full: HLos_Vector(wall, trans, recv, tunnel, simulation, &final_response, wh_coll);
    HNLos_Vector(wall, trans, recv, tunnel, simulation, &final_response, 0.2, wh_coll);
    HNLos_Vector(wall, trans, recv, tunnel, simulation, &final_response, 2.8, wh_coll);
    HScatter_Vector(wall, trans, recv, tunnel, simulation, &final_response, 2.8, wh_coll); break;
        case scattering: HScatter_Vector(wall, trans, recv, tunnel, simulation, &final_response, 2.8, wh_coll); break;
        case los: HLos_Vector(wall, trans, recv, tunnel, simulation, &final_response, wh_coll);break;
        case nlos: HNLos_Vector(wall, trans, recv, tunnel, simulation, &final_response, 0.2, wh_coll);
    HNLos_Vector(wall, trans, recv, tunnel, simulation, &final_response, 2.8, wh_coll); break;
    }

    return final_response;
}



void conv(std::vector<float> *final_response, float** h_vector, std::vector<float>* h_led, int sz){
    float sum = 0;
    for(int i=0; i<(*h_led).size(); i++){
        for(int j=0; j<sz; j++){

            (*final_response)[i+j] += (*h_led)[i]*(*h_vector)[j];
            sum += (*h_led)[i]*(*h_vector)[j];
            if((*h_vector)[j]>0 && i==0){
            }
        }

    }
}


int find_binary(std::vector<float> arr, int low, int high, float x){
    while(low<=high){
        int mid = (low+high)/2;
        if(arr[mid]<=x && arr[mid+1]>x){
            return mid;
        }else if(arr[mid]>x){
            high = mid - 1;
        }else if(arr[mid]<x){
            low = mid + 1;
        }
    }
    return -1;
}
int find_index(std::vector<float> arr, float x){
    int i = 1;
    int l = arr.size();
    while(i<l){
        if(arr[i-1]<=x && arr[i]>x){
            return i;
        }else if(2*i<l && arr[2*i-1]>x){
            return find_binary(arr,i-1,2*i,x);
        }

        i = 2*i;
    }
    return find_binary(arr, i/2, l-1, x);
}
int find_index_t(py::array_t<float> arr_t, float x){
    int i = 1;
    auto arr = std::vector<float>(arr_t.data(), arr_t.data()+arr_t.size());
    int l = arr.size();
    while(i<l){
        if(arr[i-1]<=x && arr[i]>x){
            return i;
        }else if(2*i<l && arr[2*i-1]>x){
            return find_binary(arr,i-1,2*i,x);
        }

        i = 2*i;
    }
    return find_binary(arr, i/2, l-1, x);
}
float dot_product(std::vector<float> v1, std::vector<float> v2){
    float sum = 0;
    for(int i=0; i < v1.size(); i++){
        sum += v1[i]*v2[i];
    }
    return sum;
}

float incline(float xt, float yt, float zt, float xr, float yr, float zr, float alpha, float beta, bool print = false){

    std::vector<float> path = std::vector<float>{xt-xr, yr-yt, zt-zr};
    std::vector<float> vector_normal = norm_vector_recv(alpha, 90-beta);
    std::vector<float> norm = std::vector<float>(path.size());
    std::transform(path.begin(), path.end(), norm.begin(), [](float i)->float {return i*i;});
    float sum = std::sqrt(std::reduce(norm.begin(), norm.end()));
    std::transform(path.begin(), path.end(), norm.begin(), [sum](float i)->float {return i/sum;});
    float dp = dot_product(norm, vector_normal);
  
    std::transform(vector_normal.begin(), vector_normal.end(), vector_normal.begin(), [](float i)-> float {return i*0.1;});
    if(dp<0){
        return pi;
    }

    return std::acos(dp);
}




float gain(float eta, float inct, float incr, float fov){
    float fov_radian = pi*fov/180;
    if(incr <= fov && incr>=0){
        return std::pow(eta,2)/std::sin(std::pow(fov_radian,2));
    }
    else{
        return 0;
    }
}


std::vector<float> vector_distance(float tx, float ty, float tz, float rx, float ry, float rz){
    std::vector<float> vect{rx-tx, ry-ty, rz-tz};
    return vect;
}
std::vector<float> norm_vector_trans(float alpha, float beta){
    float alphad = pi/180*alpha;
    float betad = pi/180*beta;
    std::vector<float> vect{std::cos(alphad)*std::sin(betad), std::sin(alphad)*std::sin(betad), -std::cos(betad)};
    return vect;
}
std::vector<float> norm_vector_recv(float alpha, float beta){
    float alphad = pi/180*alpha;
    float betad = pi/180*beta;
    std::vector<float> vect{std::cos(alphad)*std::sin(betad), std::sin(alphad)*std::sin(betad), std::cos(betad)};
    return vect;
}



ch::loss_and_time* HLos(float tx, float ty, float tz, float rx, float ry, float rz, float Ap, float eta, float alphat, float alphar,
    float betat, float betar, float incr, float incj, int m, float fov, float X, float Y, float t, float c, sh::WH_Probabilities* wh_coll){
    std::vector<float> dist_vector; std::vector<float> norm_vec_t; std::vector<float> norm_vec_r;
    float dist_vector_norm; float p1; float p2; float dm; float g;
    dist_vector = vector_distance(tx, ty, tz, rx, ry, rz);
    dist_vector_norm = std::sqrt(std::pow(tx-rx,2)+std::pow(ty-ry,2)+std::pow(tz-rz,2));
    norm_vec_t = norm_vector_trans(alphat, betat);
    p1 = dot_product(dist_vector, norm_vec_t);
    std::vector<float> dist_vector_neg = std::vector<float>(dist_vector.size());
    std::transform(dist_vector.begin(), dist_vector.end(), dist_vector_neg.begin(), [](float i) -> float { return -i;});
    norm_vec_r = norm_vector_recv(alphar, betar);
    p2 = dot_product(dist_vector_neg, norm_vec_r);
    g = gain(eta, incr, incj, fov);
    dm = dist_vector_norm/c;
    std::unique_ptr<sh::Shadowing_Parameters_Coll> sh_coll = std::make_unique<sh::Shadowing_Parameters_Coll>();
    sh_coll->create_collection(tx, ty, tz, rx, ry, rz, wh_coll);
    float expectancy = sh::calculate_expectancy(sh_coll.get()); 
    float poisson_mult = std::exp(-expectancy*INTENSITY*TIME);

    ch::loss_and_time* ret = new ch::loss_and_time();
    ret->loss = poisson_mult*std::abs((m+1)*Ap/(2*pi*std::pow(dist_vector_norm,2))*(std::pow(p1,m)/std::pow(dist_vector_norm,m))*(p2/dist_vector_norm)*g);

    ret->time = dm;


    return ret;
}

ch::loss_and_time* HNLos(float tx, float ty, float tz, float rx, float ry, float rz, float wx, float wy, float wz, float Aw, float pw, float Ap, float eta, float alphat, float alphar,
    float alphaw, float betat, float betar, float betaw, float incr, float incj, int m, float fov, float X, float Y, float t, float c , sh::WH_Probabilities* wh_coll){
    std::vector<float> dist_vector_tw; std::vector<float> norm_vec_t; std::vector<float> norm_vec_r;
    std::vector<float> dist_vector_wr;
    
    float dist_vector_norm_tw; float dist_vector_norm_wr; float p1; float p2; float p3; float p4; float dm; float g; std::vector<float> norm_vec_w;
    dist_vector_tw = vector_distance(wx, wy, wz, tx, ty, tz);
    std::vector<float> dist_vector_tw_neg = std::vector<float>(dist_vector_tw.size());
    std::transform(dist_vector_tw.begin(), dist_vector_tw.end(), dist_vector_tw_neg.begin(), [](float i) -> float { return -i;});
    dist_vector_norm_tw = std::sqrt(std::pow(tx-wx,2)+std::pow(ty-wy,2)+std::pow(tz-wz,2));
    dist_vector_wr = vector_distance(rx, ry, rz, wx, wy, wz);
    std::vector<float> dist_vector_wr_neg = std::vector<float>(dist_vector_wr.size());
    std::transform(dist_vector_wr.begin(), dist_vector_wr.end(), dist_vector_wr_neg.begin(), [](float i) -> float { return -i;});
    dist_vector_norm_wr = std::sqrt(std::pow(rx-wx,2)+std::pow(ry-wy,2)+std::pow(rz-wz,2));
    norm_vec_t = norm_vector_trans(alphat, betat);
    norm_vec_r = norm_vector_recv(alphar, betar);
    norm_vec_w = norm_vector_recv(90+alphaw, 90+betaw);
    p1 = dot_product(dist_vector_tw, norm_vec_t) > 0 ? dot_product(dist_vector_tw, norm_vec_t) : 0;
    p4 = dot_product(dist_vector_wr_neg, norm_vec_r);
    p2 = dot_product(dist_vector_tw_neg, norm_vec_w);
    p3 = dot_product(dist_vector_wr, norm_vec_w) > 0? dot_product(dist_vector_wr, norm_vec_w) : 0;
    g = gain(eta, incr, incj, fov);
    dm = (dist_vector_norm_wr+dist_vector_norm_tw)/c;
    std::unique_ptr<sh::Shadowing_Parameters_Coll> sh_coll_tw = std::make_unique<sh::Shadowing_Parameters_Coll>();
    sh_coll_tw->create_collection(tx, ty, tz, wx, wy, wz, wh_coll);
    std::unique_ptr<sh::Shadowing_Parameters_Coll> sh_coll_wr = std::make_unique<sh::Shadowing_Parameters_Coll>();
    sh_coll_wr->create_collection(wx, wy, wz, rx, ry, rz, wh_coll);
    float expectancy_tw = sh::calculate_expectancy(sh_coll_tw.get()); 
    float expectancy_wr = sh::calculate_expectancy(sh_coll_wr.get());
    float poisson_mult = std::exp(-(expectancy_wr+expectancy_tw)*INTENSITY*TIME);

    ch::loss_and_time* ret = new ch::loss_and_time();
    ret->loss = poisson_mult*std::abs((m+1)*Ap*Aw*pw*p1*p2*p3*p4*g/(std::pow(dist_vector_norm_tw,3+m)*
    std::pow(dist_vector_norm_wr,4)));
    ret->time = dm;

    return ret;
}

void HLos_Vector(pa::WallParameters* wall, pa::TransmitterParameters* trans, pa::ReceiverParameters* recv, pa::TunnelParameters* tunnel,
pa::SimulationParameters* simulation, std::vector<float>* final_response, sh::WH_Probabilities* wh_coll){
        float* h_vector = new float[140]();
        memset(h_vector, 0, sizeof(h_vector));

        float incidencia_tr_rad = incline(trans->coordinate[0], trans->coordinate[1], trans->coordinate[2], recv->coordinate[0], recv->coordinate[1], recv->coordinate[2], recv->alpha, recv->ele, false);
        float incidencia_tr = 180/pi*incidencia_tr_rad;
        ch::loss_and_time* hlos = HLos(trans->coordinate[0],trans->coordinate[1],trans->coordinate[2],recv->coordinate[0],recv->coordinate[1],recv->coordinate[2],
        recv->Ap, recv->eta, trans->alpha, recv->alpha, trans->beta, 90-recv->ele, incidencia_tr_rad, incidencia_tr, trans->m, recv->fov, tunnel->x, tunnel->y, simulation->t, simulation->c,
        wh_coll);
        int index = find_index(simulation->time, hlos->time);
        h_vector[index] = hlos->loss;
        conv(final_response, &h_vector, &simulation->h_led, simulation->time.size());
        free(hlos);
        delete [] h_vector;
}
void HNLos_Vector(pa::WallParameters* wall, pa::TransmitterParameters* trans, pa::ReceiverParameters* recv, pa::TunnelParameters* tunnel,
pa::SimulationParameters* simulation, std::vector<float>* final_response, float wall_pos, sh::WH_Probabilities* wh_coll){
        int a =1;
        float lx = tunnel->x; float ly = tunnel->y; float lz = tunnel->z;
        int Nx = lx*10; int Ny = ly*10; int Nz = lz*10;
        float dA = (lz*lx)/(Nx*Nz);
        float sum = 0;
        for(float kk=0; kk<lx; kk+=lx/Nx){
            for(float ll=0; ll<lz; ll+=lz/Nz){
                float r = distribution(generator);
                float s = distribution(generator);
                if(std::abs(kk-recv->coordinate[0])<=1.5 && ll>=1.5){
                    float incidencia_wr_rad = incline(kk, wall_pos, ll, recv->coordinate[0], recv->coordinate[1], recv->coordinate[2], recv->alpha, recv->ele);
                    float incidencia_wr = 180.0/pi*incidencia_wr_rad;
                    ch::loss_and_time* nhlos = HNLos(trans->coordinate[0],trans->coordinate[1],trans->coordinate[2],recv->coordinate[0],recv->coordinate[1],recv->coordinate[2],kk,wall_pos,ll
        ,dA,wall->pw,recv->Ap, recv->eta, trans->alpha, recv->alpha, r, trans->beta, 90-recv->ele,s, incidencia_wr_rad, incidencia_wr, trans->m, recv->fov, tunnel->x, tunnel->y, simulation->t,simulation->c,
        wh_coll);
                    sum += nhlos->loss;
                    float* h_vector = new float[140]();
                    memset(h_vector, 0, sizeof(h_vector));
                    int index = find_index(simulation->time, nhlos->time);
                    h_vector[index] = nhlos->loss;
                    conv(final_response, &h_vector, &simulation->h_led, simulation->time.size());
                    free(nhlos);
                    delete [] h_vector;
                }
            }
        }
}

void HScatter_Vector(pa::WallParameters* wall, pa::TransmitterParameters* trans, pa::ReceiverParameters* recv, pa::TunnelParameters* tunnel,
pa::SimulationParameters* simulation, std::vector<float>* final_response, float wall_pos, sh::WH_Probabilities* wh_coll){
    float* h_vector = new float[140]();
    memset(h_vector, 0, sizeof(h_vector));
    float N = simulation->scatters->n;
    for(int i=0; i<N; i++){
        ch::loss_and_time* hlos = HScatter(trans->coordinate[0],trans->coordinate[1],trans->coordinate[2],recv->coordinate[0],recv->coordinate[1],recv->coordinate[2],
        recv->Ap, recv->eta, trans->alpha, recv->alpha, trans->beta, 90-recv->ele, trans->m, recv->fov, tunnel->x, tunnel->y, simulation->t, simulation->c, wh_coll, simulation->scatters);
        int index = find_index(simulation->time, hlos->time);
        h_vector[index] = hlos->loss/N;
        free(hlos);
    }
    conv(final_response, &h_vector, &simulation->h_led, simulation->time.size());
    delete [] h_vector;

}
float p_mie(pa::ScatterParameters* scatters, float inc_sr){
    float g = scatters->g;
    float f = scatters->f;
    return (1-g*g/(4*pi))*std::pow((1/(1+g*g-2*g*std::cos(inc_sr))),1.5)+
    f*3*std::pow(std::cos(inc_sr),2)-0.5*std::pow((1+g*g),1.5);
}

float p_ray(pa::ScatterParameters* scatters, float inc_sr){
    float y = scatters->gamma;
    return 3*(1+3*y+(1-y)*std::pow(std::cos(inc_sr),2))/
    (16*pi*(1+2*y));

}

float gain_scattering(pa::ScatterParameters* scatters, float inc_sr){
    float pmie = p_mie(scatters, inc_sr);
    float pray = p_ray(scatters, inc_sr);
    float kr = scatters->kr;
    float km = scatters->km;
    float ks = kr + km;
    float p = scatters->p;
    float p_total = (kr/ks)*pray + (km/ks)*pmie;
    return p*p_total*std::sin(inc_sr);
}


ch::loss_and_time* HScatter(float tx, float ty, float tz, float rx, float ry, float rz, float Ap, float eta, float alphat, float alphar,
    float betat, float betar, int m, float fov, float X, float Y, float t, float c, sh::WH_Probabilities* wh_coll, pa::ScatterParameters* scatters){
        float r_sr = scatters->radius*std::rand()/RAND_MAX;
        float beta_sr = -pi/2*std::rand()/RAND_MAX + pi/2;
        float alpha_sr = -pi*std::rand()/RAND_MAX + pi;
        float xs = r_sr * std::cos(alpha_sr+alphar) * std::sin(beta_sr) + rx;
        float ys = r_sr * std::sin(alpha_sr+alphar) * std::sin(beta_sr) + ry;
        float zs = r_sr * std::cos(beta_sr) + rz;
        float inc_sr = incline(xs,ys,zs,rx,ry,rz,alphar,betar);
        std::vector<float> dist_ts = std::vector<float>{tx-xs, ty-ys, tz-zs};
        std::vector<float> dist_norm_ts_v = std::vector<float>(dist_ts.size());
        std::transform(dist_ts.begin(), dist_ts.end(), dist_norm_ts_v.begin(), [](float i)->float {return i*i;});

        float dist_norm_ts = std::sqrt(std::reduce(dist_norm_ts_v.begin(), dist_norm_ts_v.end()));
        float tot_dist = r_sr + dist_norm_ts;
        float gscattering = gain_scattering(scatters, inc_sr);
        float h_scat = (m+1)*Ap*gscattering/(2*pi*std::pow(tot_dist,2))*
        std::cos(inc_sr)*gain(eta, inc_sr, inc_sr*180/pi, fov);
        float dm = (tot_dist)/c;
        ch::loss_and_time* ret = new ch::loss_and_time();
        ret->loss = h_scat;
        ret->time = dm;
        return ret;
}

void response_subroutine(pa::WallParameters* wall, pa::TransmitterAggregate* trans, pa::ReceiverAggregate* recv, pa::TunnelParameters* tunnel,
pa::SimulationParameters* simulation, sh::WH_Probabilities* wh_coll, ResponseMode mode,  std::vector<std::vector<final_response_wrapper*>>* matrix){
    int i=0; int j = 0;

    for(auto pd: recv->receivers){
        for(auto led: trans->transmitters){
            std::unique_ptr<std::vector<float>> final_response = std::make_unique<std::vector<float>>(led_pd_response(wall, led, pd, tunnel, simulation, wh_coll, mode));
            (*matrix)[i][j] = new final_response_wrapper();

            (*matrix)[i][j]->final_response = std::move(final_response);

            j++;
        }
        j = 0;
        i++;
    }
}

void channel_subroutine(pa::WallParameters* wall, pa::TransmitterAggregate* trans, pa::ReceiverAggregate* recv, pa::TunnelParameters* tunnel,
pa::SimulationParameters* simulation, sh::WH_Probabilities* wh_coll, std::vector<std::vector<float>>* matrix){
    int i=0; int j = 0;

    for(auto pd: recv->receivers){
        for(auto led: trans->transmitters){
            (*matrix)[i][j] = led_pd_channel(wall, led, pd, tunnel, simulation, wh_coll);
            j++;
        }
        j = 0;
        i++;
    }
}

py::array_t<float> channel_calculation(pa::WallParameters* wall, pa::TransmitterAggregate* tconfig, pa::ReceiverConfigurations* rconfig,
pa::TunnelParameters* tunnel, pa::SimulationParameters* simulation, sh::WH_Probabilities* wh_coll){
    const size_t size = rconfig->receivers.size();
    std::vector<std::vector<std::vector<float>>> out_matrix(SIZE, std::vector<std::vector<float>>(4, std::vector<float>(4)));
    std::vector<std::thread> thread_vec = std::vector<std::thread>(SIZE);
    int i=0;
    for(auto recv : rconfig->receivers){
        std::thread ch(channel_subroutine,wall, tconfig, recv, tunnel, simulation, wh_coll, &out_matrix[i]);
        thread_vec[i] = std::move(ch);
        i++;
    }
    for(int k=0; k<thread_vec.size(); k++){
        thread_vec[k].join();

    }
    auto ret = py::array_t<float>({SIZE,4,4}, {4*4*sizeof(float), 4*sizeof(float), sizeof(float)});
    auto v = ret.mutable_unchecked<3>();
    for(int i=0; i<size; i++){
        for(int j=0; j<4; j++){
            for(int k=0; k<4; k++){
                v(i,j,k) = out_matrix[i][j][k];
            }
        }
    }
    return ret;

}
py::array_t<float> response_calculation_general(pa::WallParameters* wall, pa::TransmitterAggregate* tconfig, pa::ReceiverConfigurations* rconfig,
pa::TunnelParameters* tunnel, pa::SimulationParameters* simulation, sh::WH_Probabilities* wh_coll, ResponseMode mode){
    const size_t size = rconfig->receivers.size();
    std::srand(std::time(nullptr)); 

    std::cout << std::to_string(size);
    std::vector<std::vector<std::vector<final_response_wrapper*>>> out_matrix(SIZE, std::vector<std::vector<final_response_wrapper*>>(4, std::vector<final_response_wrapper*>(4)));
    std::vector<std::thread> thread_vec = std::vector<std::thread>(SIZE);
    int i=0;
    for(auto recv : rconfig->receivers){
        std::thread ch(response_subroutine, wall, tconfig, recv, tunnel, simulation, wh_coll, mode, &out_matrix[i]);
        thread_vec[i] = std::move(ch);
        
        i++;
    }
    for(int k=0; k<thread_vec.size(); k++){
        thread_vec[k].join();

    }
    auto ret = py::array_t<float>({SIZE,4,4,RESP_SIZE}, {4*4*RESP_SIZE*sizeof(float), 4*RESP_SIZE*sizeof(float), RESP_SIZE*sizeof(float), sizeof(float)});
    auto v = ret.mutable_unchecked<4>();
    for(int i=0; i<size; i++){
        for(int j=0; j<4; j++){
            for(int k=0; k<4; k++){
                std::unique_ptr<std::vector<float>> final_response_vec = std::move(out_matrix[i][j][k]->final_response); 
                std::vector<float> vec = *(final_response_vec.get());
                for(int l=0; l<RESP_SIZE; l++){
                    v(i,j,k,l) = vec[l];
                }
                delete out_matrix[i][j][k];
            }
        }
    }
    
    return ret;

}

py::array_t<float> response_calculation(pa::WallParameters* wall, pa::TransmitterAggregate* tconfig, pa::ReceiverConfigurations* rconfig,
pa::TunnelParameters* tunnel, pa::SimulationParameters* simulation, sh::WH_Probabilities* wh_coll){
    ResponseMode mode = full;
    return response_calculation_general(wall, tconfig, rconfig, tunnel, simulation, wh_coll, mode);
}
py::array_t<float> response_calculation_scattering(pa::WallParameters* wall, pa::TransmitterAggregate* tconfig, pa::ReceiverConfigurations* rconfig,
pa::TunnelParameters* tunnel, pa::SimulationParameters* simulation, sh::WH_Probabilities* wh_coll){
    ResponseMode mode = scattering;
    return response_calculation_general(wall, tconfig, rconfig, tunnel, simulation, wh_coll, mode);
}
py::array_t<float> response_calculation_los(pa::WallParameters* wall, pa::TransmitterAggregate* tconfig, pa::ReceiverConfigurations* rconfig,
pa::TunnelParameters* tunnel, pa::SimulationParameters* simulation, sh::WH_Probabilities* wh_coll){
    ResponseMode mode = los;
    return response_calculation_general(wall, tconfig, rconfig, tunnel, simulation, wh_coll, mode);
}
py::array_t<float> response_calculation_nlos(pa::WallParameters* wall, pa::TransmitterAggregate* tconfig, pa::ReceiverConfigurations* rconfig,
pa::TunnelParameters* tunnel, pa::SimulationParameters* simulation, sh::WH_Probabilities* wh_coll){
    ResponseMode mode = nlos;
    return response_calculation_general(wall, tconfig, rconfig, tunnel, simulation, wh_coll, mode);
}






