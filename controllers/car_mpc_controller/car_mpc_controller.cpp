// File:          car_mpc_controller.cpp
// Date:
// Description:
// Author:
// Modifications:

// You may need to add webots include files such as
// <webots/DistanceSensor.hpp>, <webots/Motor.hpp>, etc.
// and/or to add some other includes
#include <webots/Robot.hpp>
#include <webots/Motor.hpp>
#include <webots/Supervisor.hpp>
#include <iostream> 
#include <Eigen/Dense>
// #include "Array.hh"   //使用了MIT中的文件 其修改了向量与矩阵的名字 不然会与Eigen冲突
// #include "QuadProg++.hh"
 #include <qpOASES.hpp>
#include <fstream>

// All the webots classes are defined in the "webots" namespace
using namespace webots;
using namespace std;
using namespace Eigen;
using namespace qpOASES;
Motor *motor_FR,*motor_FL,*motor_BR,*motor_BL;
Supervisor* robot;
int testQP();
void setV(double _v,float _time_step);
//QuadProg++不好用  mit使用的是qpOASES  但是其windows下不好安装  一步切换到ubuntu

int main(int argc, char **argv) {
  //testQP();
  //Robot *robot = new Robot();
  robot = new Supervisor();   //使用Supervisor需要将机器人上supervisor设为true
  motor_FR = robot->getMotor("motor_FR");
  motor_FL = robot->getMotor("motor_FL");
  motor_BR = robot->getMotor("motor_BR");
  motor_BL = robot->getMotor("motor_BL");
  // get the time step of the current world.
  int timeStep = (int)robot->getBasicTimeStep();
  // double T = timeStep/1000;
  double T = 0.05;
  int P = 5;//预测长度

  MatrixXd Q(P,P);  //状态误差权重
  for(int i=0;i<P;i++){
    Q(i,i) = 5*1;
  }
  //std::cout << "Q  =\n" << Q  << std::endl;

  MatrixXd W(P,P);  //控制输出权重
  for(int i=0;i<P;i++){
    W(i,i) = 1;
  }

  MatrixXd Rk(P,1);  //参考值序列
  for(int i=0;i<P;i++){
    Rk(i,0) = 2;
  }

  double A_ = 1;
  double B_ = T;

  MatrixXd fei(P,1);  //AK

  MatrixXd theta(P,P);  //BK

  MatrixXd E(P,1);  //E

  MatrixXd H(P,P);  //H

  MatrixXd f(P,1);  //f
  MatrixXd g(1,P);

  MatrixXd xk(1,1);  
  MatrixXd lb(P,1);
  MatrixXd ub(P,1);

  float pos = 0;
  float vmax = 1,vmin = -1;
  
  const int num_of_elements = H.rows() * H.cols();  //元素数
  double H_matrix[num_of_elements];
  const auto k_num_of_rows = f.rows(); 
  double g_matrix[k_num_of_rows];  // NOLINT ， 与offset一般大
  const int num_of_lb = lb.rows();
  double lb_matrix[num_of_lb];
  double ub_matrix[num_of_lb];
  for(int i=0;i<num_of_lb;i++)
  {
    lb_matrix[i] = vmin;
    ub_matrix[i] = vmax;
  }

  QProblemB v_pro(P);
  Options options;
	//options.enableFlippingBounds = BT_FALSE;
	options.initialStatusBounds = ST_INACTIVE;
	options.numRefinementSteps = 1;
	options.enableCholeskyRefactorisation = 1;
	v_pro.setOptions( options );

  while (robot->step(timeStep) != -1) {
    // pos+=0.1;
    xk(0,0) = pos;
    for(int i=0 ; i<P;i++){
      fei(i,0) = pow(A_,i+1);
    }
    //  cout <<"fei "<<fei <<"\n";

    for(int i=0; i <P;++i){
      for(int j=0; j <=i;++j){
        theta(i,j)= pow(A_,i-j)*B_;
      }
    }
    //cout <<"theta"<<theta <<"\n";

    E = fei*xk-Rk;
    H = 2.0*(theta.transpose()*Q*theta+W);  //求矩阵的转置
    f = (2.0*E.transpose()*Q*theta).transpose(); 
    g = f.transpose();   //g为f的转置  
     // cout <<E <<"\n" << H <<"\n" <<f <<"\n";
     //转为qpOASES格式
    for(int i =0;i<H.rows();i++)
    {
      for(int j=0;j<H.cols();j++)
      {
        H_matrix[(i)*H.cols()+j] = H(i,j);
      }
      g_matrix[i] = g(0,i);
    }
    
    int_t nWSR = 10;
    v_pro.init(H_matrix, g_matrix, lb_matrix, ub_matrix, nWSR);
    real_t xOpt[2];
	  v_pro.getPrimalSolution( xOpt );
	  printf( "\nxOpt = [ %e, %e ];  objVal = %e\n\n", xOpt[0],xOpt[1],v_pro.getObjVal() );

    setV(xOpt[0],timeStep);
    // setV(0.5,timeStep);
    //pos = pos + xOpt[0]*T ;
    pos = robot->getSelf()->getPosition()[0];
    printf("pos:%f\n",pos);

  };

  // Enter here exit cleanup code.

  delete robot;
  return 0;
}


int testQP()
{
	/* Setup data of first QP. */
	real_t H[2 * 2] = { 1.0, 0.0, 0.0, 0.5 };
	real_t A[1 * 2] = { 1.0, 1.0 };
	real_t g[2] = { 1.5, 1.0 };
	real_t lb[2] = { 0.5, -2.0 };
	real_t ub[2] = { 5.0, 2.0 };
	real_t lbA[1] = { -1.0 };
	real_t ubA[1] = { 2.0 };
	/* Setup data of second QP. */
	real_t g_new[2] = { 1.0, 1.5 };
	real_t lb_new[2] = { 0.0, -1.0 };
	real_t ub_new[2] = { 5.0, -0.5 };
	real_t lbA_new[1] = { -2.0 };
	real_t ubA_new[1] = { 1.0 };
	/* Setting up QProblem object. */
	QProblem example(2, 1);
	/* Solve first QP. */
	int_t nWSR = 10;
	example.init(H, g, A, lb, ub, lbA, ubA, nWSR);

	real_t xOpt[2];
	real_t yOpt[2+1];
	example.getPrimalSolution( xOpt );
	example.getDualSolution( yOpt );
	printf( "\nxOpt = [ %e, %e ];  yOpt = [ %e, %e, %e ];  objVal = %e\n\n", 
			xOpt[0],xOpt[1],yOpt[0],yOpt[1],yOpt[2],example.getObjVal() );

	/* Solve second QP. */
	nWSR = 10;
	example.hotstart(g_new, lb_new, ub_new, lbA_new, ubA_new, nWSR);
	/* Get and print solution of second QP. */
//	real_t xOpt[2];
	example.getPrimalSolution(xOpt);
	printf("\n xOpt = [ %e, %e ]; objVal = %e\n\n",
		xOpt[0], xOpt[1], example.getObjVal());
	return 0;
}

float wheel_p=0;
void setV(double _v,float _time_step)
{
  double r = 0.05;
  _v = -_v;
  wheel_p += _time_step/1000*_v/r;
  motor_FR->setPosition(wheel_p);
  motor_FL->setPosition(wheel_p);  
  motor_BR->setPosition(wheel_p);
  motor_BL->setPosition(wheel_p);  
  

}