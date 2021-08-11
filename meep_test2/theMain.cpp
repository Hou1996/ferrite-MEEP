#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <math.h>

#include <meep.hpp>

using namespace meep;

double eps(const vec &v) { return 1.0; }

int main(int argc, char **argv)
{
  std::ofstream mon_po("monitor.dat"); //存储观察点位置的数据
                                       //  std::ofstream source_ou("sources.dat");
  initialize mpi(argc, argv);          // do this even for non-MPI Meep

  double nx = 2.3, ny = 1.02, nz = 13.0, unit_num = 20, fr = 0.31, pml_l = 2.0;
  /*std::cout << "模型以原点为坐标,单位长度1mm" << std::endl;
  std::cout << "输入模型x方向长度" << std::endl;
  std::cin >> nx;
  std::cout << "输入模型y方向长度" << std::endl;
  std::cin >> ny;
  std::cout << "输入模型z方向长度" << std::endl;
  std::cin >> nz;
  std::cout << "输入模型单位长度像素数量" << std::endl;
  std::cin >> unit_num;
  std::cout << "输入频率" << std::endl;
  std::cin >> fr;
  std::cout << "输入PML厚度" << std::endl;
  std::cin >> pml_l;
  double match_cw, match_cl, plus_cw, plus_cl, fer_w, fer_l, fer_aw, waveguide_Y;
  std::cout << "输入介质区域的模型参数" << std::endl;
  std::cout << "输入匹配陶瓷宽度" << std::endl;
  std::cin >> match_cw;
  std::cout << "输入匹配陶瓷长度" << std::endl;
  std::cin >> match_cl;
  std::cout << "输入加载陶瓷宽度" << std::endl;
  std::cin >> plus_cw;
  std::cout << "输入加载陶瓷长度" << std::endl;
  std::cin >> plus_cl;
  std::cout << "输入铁氧体宽度" << std::endl;
  std::cin >> fer_w;
  std::cout << "输入铁氧体长度" << std::endl;
  std::cin >> fer_l;
  std::cout << "输入铁氧体磁环中间空气孔宽度" << std::endl;
  std::cin >> fer_aw;
  std::cout << "波导的高(Y)" << std::endl;
  std::cin >> waveguide_Y;*/

  double resolution = unit_num;                  // pixels per distance
  grid_volume v = vol3d(nx, ny, nz, resolution); // 三维模型:nx * ny * nz 3d cell
  material_function s2;
  structure s(v, eps, pml(pml_l, Z)); // 计算区域各区域材料参数设置
  vec cen = v.center();

  std::cout << "set material" << std::endl;
  double feps = 14.0, leps = 35.0, meps = 6.5;
  //double feps, leps, meps;
  /*  std::cout << "输入铁氧体相对介电常数" << std::endl;
  std::cin >> feps;
  std::cout << "输入加载陶瓷相对介电常数 " << std::endl;
  std::cin >> leps;
  std::cout << "输入匹配陶瓷相对介电常数" << std::endl;
  std::cin >> meps;*/
  double decay = 0.02, remanence_r = 0.75, remanence = 4000, satmag = 1.12E5;
  /*  double decay, remanence_r, remanence, satmag;
  std::cout << "输入衰减系数" << std::endl;
  std::cin >> decay;
  std::cout << "输入剩磁比" << std::endl;
  std::cin >> remanence_r;
  std::cout << "输入剩磁化强度" << std::endl;
  std::cin >> remanence;
  std::cout << "输入饱和磁化强度" << std::endl;
  std::cin >> satmag;*/
  double free = 9.3E9;
  /* double free;
  std::cout << "输入材料中应用的频率(真实物理空间频率)" << std::endl;
  std::cin >> free;*/
  ferrite_material_function fm; //
  fm.set_center(cen);
  fm.set_freq(free);
  fm.set_model(0.3, 0.5, 0.2, 6., 0.1, 6., 0.2, 1.02);
  //  fm.set_model(match_cw, match_cl, plus_cw, plus_cl, fer_w, fer_l, fer_aw, waveguide_Y);
  fm.set_eps(feps, leps, meps);
  fm.set_mu_pare(decay, remanence_r, remanence, satmag);

  //设置磁导率张量矩阵，参数为空：计算模型的张量矩阵， 参数为0：计算磁导率矩阵为单位矩阵
  int a = 2;
  fm.set_muMat(0);

  s.set_materials(fm, false);
  //用于设置输入平面，输出平面的几个坐标点
  vec p1(0, 0, pml_l + 0.2);
  vec p2(nx, ny, pml_l + 0.2);
  vec p3(0, 0, nz / 2);
  vec p4(nx, ny, nz / 2);
  vec p5(0, ny / 2, 0);
  vec p6(nx, ny / 2, nz);
  vec p7(0, 0, nz - pml_l - 0.2);
  vec p8(nx, ny, nz - pml_l - 0.2);

  //field set
  fields f(&s);
  //field boundary
  f.set_boundary(High, X, Metallic);
  f.set_boundary(Low, X, Metallic);
  f.set_boundary(High, Y, Metallic);
  f.set_boundary(Low, Y, Metallic);

  // output real(eps) / real(mu)
  f.output_hdf5(Dielectric, volume(p5, p6));
  f.output_hdf5(Permeability, volume(p5, p6));

  // source f = f_0 / c_0
  double powerof = 5.E3;
  //std::cin>>powerof;
  gaussian_src_time src(fr, 0.01);                  //gauss*sin  参数：中心频率， 频宽  阅读源码有其他重载函数
  continuous_src_time csrc(fr, 0.32, 0, 100);       //(cos,sin)
  square_pulse_time sq_src{powerof, fr, 100, 0, 0}; // poweof * sin * square_pulse
  std::complex<double> amp = {1.0, 1.0};            //振幅
  //f.add_point_source(Ex, src2, vec(nx / 2, ny / 2, nz / 2), amp);  //点源

  //  f.add_volume_source(Ey, csrc, volume(p3, p4), amp);
  f.add_volume_source(Ey, src, volume(p1, p2), amp); //极化方向为Ey，波型为高斯调制函数，所在平面是(平面左下角坐标，平面右上角坐标)，振幅是amp
  f.add_volume_source(Hz, src, volume(p1, p2), amp); //极化方向为Hz，波型为高斯调制函数，所在平面是(平面左下角坐标，平面右上角坐标)，振幅是amp

  //
  monitor_point m1, m2;

  std::cout << f.last_source_time() << std::endl;
  int n = 0;
  // f.t = 6000;
  while (f.time() < f.last_source_time() * 20)
  {
    //    source_ou << f.time() << "  " << std::real(src.dipole(f.time())) << "  " << std::imag(src.dipole(f.time())) << std::endl;
    f.step();

    //设置观察点坐标
    f.get_point(&m1, vec(nx / 2, ny / 2, pml_l + 1));
    f.get_point(&m2, vec(nx / 2, ny / 2, nz - pml_l - 1));
    //输出观察点位置电场各方向分量
    mon_po << f.time()
           << "     " << std::real(m1.get_component(Ex)) << "  " << std::imag(m1.get_component(Ex))
           << "     " << std::real(m1.get_component(Ey)) << "  " << std::imag(m1.get_component(Ey))
           << "     " << std::real(m1.get_component(Ez)) << "  " << std::imag(m1.get_component(Ez))
           << "     " << std::real(m2.get_component(Ex)) << "  " << std::imag(m2.get_component(Ex))
           << "     " << std::real(m2.get_component(Ey)) << "  " << std::imag(m2.get_component(Ey))
           << "     " << std::real(m2.get_component(Ez)) << "  " << std::imag(m2.get_component(Ez))
           << std::endl;

    //输出特定时刻的平面volume(p5,p6)de 电磁场量场分布
    //    if (n > 5000)
    {
      if ((n % 4000) == 0)
        f.output_hdf5(Ey, volume(p5, p6));
      if ((n % 4000) == 0)
        f.output_hdf5(Hx, volume(p5, p6));
      if ((n % 4000) == 0)
        f.output_hdf5(Hz, volume(p5, p6));
    }
    n++;
  }

  //计算结束输出平面电磁场量分布图
  f.output_hdf5(Ex, volume(p5, p6));
  f.output_hdf5(Ey, volume(p5, p6));
  f.output_hdf5(Ez, volume(p5, p6));
  f.output_hdf5(Hx, volume(p5, p6));
  f.output_hdf5(Hy, volume(p5, p6));
  f.output_hdf5(Hz, volume(p5, p6));

  return 0;
}