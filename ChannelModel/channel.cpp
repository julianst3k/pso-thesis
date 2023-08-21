#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <vector>
#include <cmath>
#include <algorithm>
#include <parameter_aggregate.h>
#include <random>
#include <channel.h>

#define pi 3.1415926535
namespace pa = parameter_aggregate;
namespace py = pybind11;

void initialize_3d_matrix(float*** out_matrix, int size, int h, int w){
    out_matrix = new float**[size];
    for(int x=0; x<size; x++){
        out_matrix[x] = new float*[h];
        for(int y=0; y<h; y++){
            out_matrix[x][y] = new float[w];
        }
    }
}

float led_pd_channel(pa::WallParameters* wall, pa::TransmitterParameters* trans, pa::ReceiverParameters* recv, pa::TunnelParameters* tunnel,
pa::SimulationParameters* simulation){
    std::vector<float> final_response(simulation->time.size()+simulation->h_led.size()-1);
    HLos_Vector(wall, trans, recv, tunnel, simulation, final_response);
    HNLos_Vector(wall, trans, recv, tunnel, simulation, final_response);
    HNLos_Vector(wall, trans, recv, tunnel, simulation, final_response);
    float sum = std::reduce(std::cbegin(final_response), std::cend(final_response));
    final_response.clear();
    return sum;
}



void conv(std::vector<float> final_response, float* h_vector, std::vector<float> h_led, int sz){
    for(int i=0; i<h_led.size(); i++){
        for(int j=0; j<sz; j++){
            final_response[i+j] = h_led[i]*h_vector[j];
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
    auto arr = std::vector<T>(arr_t.data(), arr_t.data()+arr_t.size());
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
float incline(float xt, float yt, float zt, float xr, float yr, float zr, float alpha, float beta){
    std::vector<float> path = std::vector<float>{xr-xt, yr-yt, zt-zr};
    std::vector<float> vector_normal = norm_vector_recv(alpha, 90-beta);
    std::vector<float> norm;
    std::transform(path.begin(), path.end(), norm.begin(), [](float i)->float {return i*i;});
    float sum = std::reduce(norm.begin(), norm.end());
    std::transform(norm.begin(), norm.end(), norm.begin(), [sum](float i)->float {return i/sum;});
    float dp = dot_product(norm, vector_normal);
    std::transform(vector_normal.begin(), vector_normal.end(), vector_normal.begin(), [](float i)-> float {return i*0.1;});
    if(dp>0){
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
    std::vector<float> vect{-std::cos(alphad)*std::sin(betad), -std::sin(alphad)*std::sin(betad), std::cos(betad)};
    return vect;
}



ch::loss_and_time* HLos(float tx, float ty, float tz, float rx, float ry, float rz, float Ap, float eta, float alphat, float alphar,
    float betat, float betar, float incr, float incj, int m, float fov, float X, float Y, float t, float c){
    std::vector<float> dist_vector; std::vector<float> norm_vec_t; std::vector<float> norm_vec_r;
    std::vector<float> dist_vector_neg;
    float dist_vector_norm; float p1; float p2; float dm; float g;
    dist_vector = vector_distance(tx, ty, tz, rx, ry, rz);
    dist_vector_norm = std::sqrt(std::pow(tx-rx,2)+std::pow(ty-ry,2)+std::pow(tz-rz,2));
    norm_vec_t = norm_vector_trans(alphat, betat);
    p1 = dot_product(dist_vector, norm_vec_t);
    std::transform(dist_vector.begin(), dist_vector.end(), dist_vector_neg.begin(), [](float i) -> float { return -i;});
    norm_vec_r = norm_vector_recv(alphar, betar);
    p2 = dot_product(dist_vector_neg, norm_vec_r);
    g = gain(eta, incr, incj, fov);
    dm = dist_vector_norm/c;
    ch::loss_and_time* ret = new ch::loss_and_time();
    ret->loss = std::abs((m+1)*Ap/(2*pi*std::pow(dist_vector_norm,2))*(std::pow(p1,m)/dist_vector_norm)*(p2/dist_vector_norm)*g);
    ret->time = dm;

    return ret;
}

ch::loss_and_time* HNLos(float tx, float ty, float tz, float rx, float ry, float rz, float wx, float wy, float wz, float Aw, float pw, float Ap, float eta, float alphat, float alphar,
    float alphaw, float betat, float betar, float betaw, float incr, float incj, int m, float fov, float X, float Y, float t, float c){
    std::vector<float> dist_vector_tw; std::vector<float> norm_vec_t; std::vector<float> norm_vec_r;
    std::vector<float> dist_vector_tw_neg; std::vector<float> dist_vector_wr; std::vector<float> dist_vector_wr_neg;
    float dist_vector_norm_tw; float dist_vector_norm_wr; float p1; float p2; float p3; float p4; float dm; float g; std::vector<float> norm_vec_w;
    dist_vector_tw = vector_distance(tx, ty, tz, wx, wy, wz);
    std::transform(dist_vector_tw.begin(), dist_vector_tw.end(), dist_vector_tw_neg.begin(), [](float i) -> float { return -i;});
    dist_vector_norm_tw = std::sqrt(std::pow(tx-wx,2)+std::pow(ty-wy,2)+std::pow(tz-wz,2));
    dist_vector_wr = vector_distance(wx, wy, wz, rx, ry, rz);
    std::transform(dist_vector_wr.begin(), dist_vector_wr.end(), dist_vector_wr_neg.begin(), [](float i) -> float { return -i;});
    dist_vector_norm_wr = std::sqrt(std::pow(rx-wx,2)+std::pow(ry-wy,2)+std::pow(rz-wz,2));
    norm_vec_t = norm_vector_trans(alphat, betat);
    norm_vec_r = norm_vector_recv(alphar, betar);
    norm_vec_w = norm_vector_recv(alphaw, betaw);
    p1 = dot_product(dist_vector_tw, norm_vec_t);
    p4 = dot_product(dist_vector_wr_neg, norm_vec_r);
    p2 = dot_product(dist_vector_tw_neg, norm_vec_w);
    p3 = dot_product(dist_vector_wr, norm_vec_w);
    g = gain(eta, incr, incj, fov);
    dm = (dist_vector_norm_wr+dist_vector_norm_tw)/c;

    ch::loss_and_time* ret = new ch::loss_and_time();
    ret->loss = std::abs((m+1)*Ap*Aw*pw*p1*p2*p3*p4*g/(std::pow(dist_vector_norm_tw,4)*
    std::pow(dist_vector_norm_wr,4)*pi*pi));
    ret->time = dm;

    return ret;
}

void HLos_Vector(pa::WallParameters* wall, pa::TransmitterParameters* trans, pa::ReceiverParameters* recv, pa::TunnelParameters* tunnel,
pa::SimulationParameters* simulation, std::vector<float> final_response){
        float* h_vector;
        h_vector = new float[simulation->time.size()];
        float incidencia_tr = incline(trans->coordinate[0], trans->coordinate[1], trans->coordinate[2], recv->coordinate[0], recv->coordinate[1], recv->coordinate[2], recv->alpha, recv->ele);
        float incidencia_tr_rad = pi/180*incidencia_tr;
        ch::loss_and_time* hlos = HLos(trans->coordinate[0],trans->coordinate[1],trans->coordinate[2],recv->coordinate[0],recv->coordinate[1],recv->coordinate[2],
        recv->Ap, recv->eta, trans->alpha, recv->alpha, trans->beta, recv->ele, incidencia_tr_rad, incidencia_tr, trans->m, recv->fov, tunnel->x, tunnel->y, simulation->t, simulation->c);
        int index = find_index(simulation->time, hlos->time);
        h_vector[index] = hlos->loss;
        conv(final_response, h_vector, simulation->h_led, simulation->time.size());
        free(h_vector);
}
void HNLos_Vector(pa::WallParameters* wall, pa::TransmitterParameters* trans, pa::ReceiverParameters* recv, pa::TunnelParameters* tunnel,
pa::SimulationParameters* simulation, std::vector<float> final_response){
        int a =1;
        float lx = tunnel->x; float ly = tunnel->y; float lz = tunnel->z;
        int Nx = lx*3; int Ny = ly*3; int Nz = lz*3;
        float dA = (lz*lx)/(Nx*Nz);
        for(float kk=0; kk<lx; kk+=lx/Nx){
            for(float ll=0; ll<ly; ll+=ly/Ny){
                std::vector<float> vec = std::vector<float>{kk, 0.2, ll};
                std::default_random_engine generator;
                std::uniform_real_distribution<float> distribution(0.0,90);
                float r = distribution(generator);
                float s = distribution(generator);
                if(std::abs(kk-recv->coordinate[0])<=1.5 && ll>=1.5){
                    float incidencia_wr = incline(kk, 0.2, ll, recv->coordinate[0], recv->coordinate[1], recv->coordinate[2], recv->alpha, recv->ele);
                    float incidencia_wr_rad = pi/180*incidencia_wr;
                    ch::loss_and_time* nhlos = HNLos(trans->coordinate[0],trans->coordinate[1],trans->coordinate[2],kk,0.2,ll,recv->coordinate[0],recv->coordinate[1],recv->coordinate[2]
        ,dA/70,wall->pw,recv->Ap, recv->eta, trans->alpha, r, recv->alpha, trans->beta, s, recv->ele, incidencia_wr_rad, incidencia_wr, trans->m, recv->fov, tunnel->x, tunnel->y, simulation->t,simulation->c);
                    float* h_vector;
                    h_vector = new float[simulation->time.size()];
                    int index = find_index(simulation->time, nhlos->time);
                    h_vector[index] = nhlos->loss;
                    conv(final_response, h_vector, simulation->h_led, simulation->time.size());
                    free(h_vector);
                }
            }
        }
}

void channel_subroutine(pa::WallParameters* wall, pa::TransmitterAggregate* trans, pa::ReceiverAggregate* recv, pa::TunnelParameters* tunnel,
pa::SimulationParameters* simulation, std::vector<std::vector<float>> matrix){
    int i=0; int j = 0;

    for(auto pd: recv->receivers){
        for(auto led: trans->transmitters){
            matrix[i][j] = led_pd_channel(wall, led, pd, tunnel, simulation);
            j++;
        }
        i++;
    }
}

py::array_t<float> channel_calculation(pa::WallParameters* wall, pa::TransmitterAggregate* tconfig, pa::ReceiverConfigurations* rconfig,
pa::TunnelParameters* tunnel, pa::SimulationParameters* simulation){
    const size_t size = rconfig->receivers.size();
    std::vector<std::vector<std::vector<float>>> out_matrix(16, std::vector<std::vector<float>>(4, std::vector<float>(4)));
    int i=0;
    for(auto recv : rconfig->receivers){
        channel_subroutine(wall, tconfig, recv, tunnel, simulation, out_matrix[i]);
        i++;
    }
    float* out_matrix_flat = (float*)malloc(16*16*sizeof(float));
    int cnt = 0;
    for(int i=0; i<size; i++){
        for(int j=0; j<4; j++){
            for(int k=0; k<4; k++){
                out_matrix_flat[cnt] = out_matrix[i][j][k];
                cnt++;
            }
        }
    }
    return py::array_t<float>({16,4,4}, {4*4*sizeof(float), 4*sizeof(float), sizeof(float)}, out_matrix_flat);

}



