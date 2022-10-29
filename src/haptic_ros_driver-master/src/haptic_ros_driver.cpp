#include "ros/ros.h"
#include "std_msgs/String.h"
#include <geometry_msgs/Wrench.h>

#include "haptic_ros_driver/HapticDevice.h"
#include <ur_rtde/rtde_control_interface.h>
#include <ur_rtde/rtde_receive_interface.h>
#include <thread>
#include <chrono>
#include "haptic_ros_driver/dhdc.h"
#include <Eigen/Dense>
#include <cmath>
#include "iir_filters/Iir.h"
#include <fstream>

// new include
#include <Eigen/Dense>
#include <atomic>
#include <chrono>
#include <kdl/chainiksolvervel_pinv.hpp>
#include <kdl/frames.hpp>
#include <kdl/jntarray.hpp>
#include <plog/Appenders/ColorConsoleAppender.h>
#include <plog/Initializers/RollingFileInitializer.h>
#include <plog/Log.h>
#include <ros/ros.h>
#include <sensor_msgs/JointState.h>
#include <std_msgs/Int64.h>
#include <std_msgs/Int64MultiArray.h>
#include <thread>
#include <trac_ik/trac_ik.hpp>
#include <ur_rtde/rtde_control_interface.h>
#include <ur_rtde/rtde_io_interface.h>
#include <ur_rtde/rtde_receive_interface.h>
#include <velocity_profile/interpolate.h>

using namespace ur_rtde;
using namespace std::chrono;
using namespace KDL;
using namespace std;

/*active compliance*/
/**
 * @brief 四元素插值公式
 *
 * @param end 终止四元素
 * @param s 百分比
 * @return std::vector< double > 在s位置处的四元素
 */
std::vector<double> UnitQuaternion_intep(const std::vector<double> &start, const std::vector<double> &end, double s)
{
	if (s > 1 || s < 0)
	{
		std::cerr << "values of S outside interval [0,1]" << std::endl;
	}

	double cosTheta = start[0] * end[0] + start[1] * end[1] + start[2] * end[2] + start[3] * end[3];
	std::vector<double> start_2 = start;

	//** 这里是为了取最短路径 **//
	if (cosTheta < 0)
	{
		for (int i = 0; i < 4; i++)
			start_2[i] *= -1;
		cosTheta *= -1;
	}
	//**-------------------------------**//

	double theta = acos(cosTheta);
	if (theta == 0 || s == 0)
		return start_2;
	else
	{
		double coefficient_1 = sin((1 - s) * theta) / sin(theta);
		double coefficient_2 = sin((s)*theta) / sin(theta);

		return std::vector<double>{
			coefficient_1 * start_2[0] + coefficient_2 * end[0],
			coefficient_1 * start_2[1] + coefficient_2 * end[1],
			coefficient_1 * start_2[2] + coefficient_2 * end[2],
			coefficient_1 * start_2[3] + coefficient_2 * end[3]};
	}
}

/**
 * @brief 直线规划（使用doubleS速度曲线）
 *
 * @param start 起始位姿
 * @param end 终止位姿
 * @return std::vector< KDL::Frame > 轨迹
 */
std::vector<KDL::Frame> moveL(KDL::Frame start, KDL::Frame end, double path_lenth = 0)
{
	//** 变量初始化 **//
	KDL::Vector Pstart = start.p;
	KDL::Vector Pend = end.p;
	std::vector<KDL::Frame> traj_1;
	double s = 0;
	rocos::R_INTERP doubleS;
	bool isplanned = doubleS.planDoubleSPorfile(0, 0, 1, 0, 0, 0.1 / path_lenth, 1 / path_lenth, 6 / path_lenth);
	const double timegap = 0.002; // 500hz
	assert(isplanned && (doubleS.getDuration() > 0) && "doubleS planning is failed");
	int N = doubleS.getDuration() / timegap;
	std::vector<double> Quaternion_start{0, 0, 0, 0};
	std::vector<double> Quaternion_end{0, 0, 0, 0};
	std::vector<double> Quaternion_interp{0, 0, 0, 0};
	start.M.GetQuaternion(Quaternion_start.at(0), Quaternion_start.at(1), Quaternion_start.at(2), Quaternion_start.at(3));
	end.M.GetQuaternion(Quaternion_end.at(0), Quaternion_end.at(1), Quaternion_end.at(2), Quaternion_end.at(3));

	//**-------------------------------**//

	//** 轨迹计算 **//
	for (int i = 0; i <= N; i++)
	{
		s = doubleS.pos(timegap * i);
		KDL::Vector P = Pstart + (Pend - Pstart) * s;
		Quaternion_interp = UnitQuaternion_intep(Quaternion_start, Quaternion_end, s);
		traj_1.push_back(KDL::Frame(KDL::Rotation::Quaternion(Quaternion_interp[0], Quaternion_interp[1], Quaternion_interp[2], Quaternion_interp[3]), P));
	}
	//**-------------------------------**//
	return traj_1;
}

/**
 * @brief 计算直线轨迹的速度向量
 *
 * @param Cartesian_vel 计算的笛卡尔速度向量
 * @param t 求解时间点
 * @param start 起始位姿
 * @param end 终止位姿
 * @return int
 */
int moveL_vel(KDL::Twist &Cartesian_vel, double t, KDL::Frame start, KDL::Frame end)
{
	KDL::Vector vel = end.p - start.p;
	double path_lenth = vel.Normalize();

	rocos::R_INTERP doubleS;
	bool isplanned = doubleS.planDoubleSPorfile(0, 0, 1, 0, 0, 0.1 / path_lenth, 1 / path_lenth, 6 / path_lenth);
	vel = vel * doubleS.vel(t);
	Cartesian_vel.vel = vel;

	KDL::Rotation rot_start_end = start.M.Inverse() * end.M;
	KDL::Vector aixs = rot_start_end.GetRot();
	aixs.Normalize();
	Cartesian_vel.rot = aixs * doubleS.vel(t);

	return 0;
}

class ur_force_control
{
public:
	static const unsigned int joint_num{6};
	//机器人读取变量
private:
	KDL::JntArray max_joints_vel{joint_num};

	std::vector<double> TCP_force{0, 0, 0};
	std::vector<double> force_pos_offset{0, 0, 0};
	std::vector<double> force_vel_offset{0, 0, 0};
	std::vector<double> force_last_vel_offset{0, 0, 0};
	std::vector<double> force_acc_offset{0, 0, 0};
	std::vector<double> force_last_acc_offset{0, 0, 0};

	std::vector<double> TCP_torque{0, 0, 0};
	std::vector<double> torque_pos_offset{0, 0, 0};
	std::vector<double> torque_vel_offset{0, 0, 0};
	std::vector<double> torque_last_vel_offset{0, 0, 0};
	std::vector<double> torque_acc_offset{0, 0, 0};
	std::vector<double> torque_last_acc_offset{0, 0, 0};

	std::vector<double> M{30, 30, 30, 30, 30, 30};
	std::vector<double> K{500, 500, 500, 500, 500, 500};
	std::vector<double> B{30, 30, 30, 30, 30, 30};

	KDL::Twist _Cartesian_vel;
	double dt{0.002};

	std::shared_ptr<RTDEControlInterface> _ur_control_ptr{nullptr};

public:
	ur_force_control(std::shared_ptr<RTDEControlInterface> ur_control_ptr)
		: _ur_control_ptr(ur_control_ptr)
	{
		for (int i{0}; i < joint_num; i++)
			B[i] = 2 * 1 * sqrt(M[i] * K[i]);

		for (int i{0}; i < joint_num; i++)
			max_joints_vel(i) = 40 * M_PI / 180;
	}

	void calculate_translate()
	{
		for (int i{0}; i < 3; i++)
		{
			force_acc_offset[i] = (TCP_force[i] - B[i] * force_vel_offset[i] - K[i] * force_pos_offset[i]) / M[i];
			force_vel_offset[i] = dt * (force_acc_offset[i] + force_last_acc_offset[i]) / 2 + force_vel_offset[i];
			force_pos_offset[i] = dt * (force_vel_offset[i] + force_last_vel_offset[i]) / 2 + force_pos_offset[i];

			force_last_acc_offset[i] = force_acc_offset[i];
			force_last_vel_offset[i] = force_vel_offset[i];

			_Cartesian_vel.vel[i] = force_vel_offset[i];
		}
	}

	KDL::Rotation calculate_rotation()
	{
		KDL::Vector delta_rot;
		KDL::Vector current_rot;
		KDL::Rotation template_rot;

		for (int i{0}; i < 3; i++)
		{
			torque_acc_offset[i] = (TCP_torque[i] - B[i] * torque_vel_offset[i] - K[i] * torque_pos_offset[i]) / M[i];
			torque_vel_offset[i] = dt * (torque_acc_offset[i] + torque_last_acc_offset[i]) / 2 + torque_vel_offset[i];

			delta_rot(i) = dt * (torque_vel_offset[i] + torque_last_vel_offset[i]) / 2;
			current_rot(i) = torque_pos_offset[i];

			torque_last_acc_offset[i] = torque_acc_offset[i];
			torque_last_vel_offset[i] = torque_vel_offset[i];

			_Cartesian_vel.rot[i] = torque_vel_offset[i];
		}

		template_rot = KDL::Rotation::Rot(delta_rot, delta_rot.Norm()) * KDL::Rotation::Rot(current_rot, current_rot.Norm());

		current_rot = template_rot.GetRot();

		for (int i{0}; i < 3; i++)
			torque_pos_offset[i] = current_rot[i];

		return template_rot;
	}

	int calculate(KDL::Frame &pos_offset, KDL::Twist &Cartesian_vel)
	{
		calculate_translate();

		for (int i{0}; i < 3; i++)
		{
			pos_offset.p[i] = force_pos_offset[i];
		}

		pos_offset.M = calculate_rotation();

		Cartesian_vel = _Cartesian_vel;

		return 0;
	}

	void check_before_move(KDL::JntArray &joints_vel)
	{
		for (int i{0}; i < joint_num; i++)
		{
			if (joints_vel(i) > max_joints_vel(i) || TCP_force[i > 2 ? 2 : i] > 70 || TCP_torque[i > 2 ? 2 : i] > 30)
			{
				PLOG_ERROR << "joint [" << i << "]  速度过快";
				PLOG_ERROR << "目标速度  = " << joints_vel(i);
				PLOG_ERROR << "允许最大速度  = " << max_joints_vel(i);

				if (_ur_control_ptr)
				{
					_ur_control_ptr->servoStop();
					_ur_control_ptr->stopScript();
				}
				exit(0);
			}
		}
	}

	void set_force(double force_x, double force_y, double force_z)
	{
		TCP_force[0] = force_x;
		TCP_force[1] = force_y;
		TCP_force[2] = force_z;
		// PLOG_DEBUG.printf( "TCP_force  = %f %f %f", TCP_force[ 0 ], TCP_force[ 1 ], TCP_force[ 2 ] );
	}

	void set_torque(double tor_que_x, double tor_que_y, double tor_que_z)
	{
		TCP_torque[0] = tor_que_x;
		TCP_torque[1] = tor_que_y;
		TCP_torque[2] = tor_que_z;
		// PLOG_DEBUG.printf( "TCP_torque  = %f %f %f", TCP_torque[ 0 ], TCP_torque[ 1 ], TCP_torque[ 2 ] );
	}

	void set_damp(double value)
	{
		static double damp = 1;
		damp += value;

		for (int i{0}; i < joint_num; i++)
			B[i] = 2 * damp * sqrt(M[i] * K[i]);

		PLOG_DEBUG << "damp  = " << damp;
	}
};

void RunTeleop_damp(const std_msgs::Int64::ConstPtr &msg, ur_force_control *ur_admittance_ptr)
{
	if (msg->data == 9)
		ur_admittance_ptr->set_damp(0.1);
	else if (msg->data == 10)
		ur_admittance_ptr->set_damp(-0.1);
}

void RunTeleop_force_torque(const std_msgs::Int64MultiArray::ConstPtr &msg, ur_force_control *ur_admittance_ptr)
{
	ur_admittance_ptr->set_force(msg->data[0], msg->data[1], msg->data[2]);
	ur_admittance_ptr->set_torque(msg->data[3], msg->data[4], msg->data[5]);
}

/*hhhh*/

/**
 * This tutorial demonstrates simple sending of messages over the ROS system.
 */
int main(int argc, char *argv[])
{

	ros::init(argc, argv, "haptic_ros_driver");

	ros::NodeHandle nh("~");

	//连接UR
	RTDEControlInterface rtde_control("192.168.3.101");
	RTDEReceiveInterface rtde_receive("192.168.3.101");

	HapticDevice haptic_dev(nh, false); //实例化HapticDevice，haptic_dev是一个HapticDevice类的对象
	/**创建Publisher**/
	/**********/
	std::vector<double> TcpPose = rtde_receive.getActualTCPPose();
	std::vector<double> TcpPoseInit = rtde_receive.getActualTCPPose();

	int button0_state_ = dhdGetButton(0);
	double PositonScale = 1.0; //位置的scale
	double PoseScale = 1.0;	   //姿态的scale
	haptic_dev.Start();

	ros::AsyncSpinner spinner(3);
	spinner.start();

	haptic_dev.keep_alive_ = true;

	double current_position[3] = {0.0, 0.0, 0.0};
	double current_position_Init[3] = {0.0, 0.0, 0.0};
	double current_orientRad[3] = {0.0, 0.0, 0.0};
	double current_orientRad_Init[3] = {0.0, 0.0, 0.0};
	double Omega_matrix[3][3];
	double Omega_matrixInit[3][3];
	Eigen::Matrix3d Omega_matrixEigen;
	Eigen::Matrix3d Omega_matrixInitEigen;
	Eigen::Matrix3d DeltaR;
	Eigen::Matrix3d Tcp_Rotation;
	Eigen::Vector3d axis;
	double alpha{0.0};
	Eigen::Vector3d axisInit;
	double alphaInit{0.0};
	Eigen::Matrix3d Tcp_Rotation_Init;
	Tcp_Rotation_Init.setIdentity();
	Eigen::Matrix3d R;
	Eigen::AngleAxisd K;
	Eigen::Vector3d Khat;
	std::vector<double> TcpForce(6, 0.0);
	rtde_control.zeroFtSensor(); //将力传感器的读数设为0
	TcpForce = rtde_receive.getActualTCPForce();

	std::cout << TcpForce[0] << " " << TcpForce[1] << " " << TcpForce[2]
			  << " " << TcpForce[1] << " " << TcpForce[4] << " " << TcpForce[5] << std::endl;
	// filter function
	//  init filter
	float fx, fy, fz, tx, ty, tz;
	float scaling = 1;									 //力的scale
	Iir::Butterworth::LowPass<2> f1, f2, f3, f4, f5, f6; // NOTE： here order should replaced by a int number!
	const float samplingrate = 200;						 // Hz
	const float cutoff_frequency = 0.5;					 // Hz
	f1.setup(2, samplingrate, cutoff_frequency);		 // NOTE： here order should replaced by a int number!
	f2.setup(2, samplingrate, cutoff_frequency);		 // NOTE： here order should replaced by a int number!
	f3.setup(2, samplingrate, cutoff_frequency);		 // NOTE： here order should replaced by a int number!

	// Set input parameters

	ur_force_control ur_admittance{nullptr};
	KDL::Frame frame_offset;
	KDL::Twist admittance_vel;

	std::vector<double> std_tcp_pos_command(6, 0);

#if 1

	while (ros::ok() && (haptic_dev.keep_alive_ == true))
	{
		button0_state_ = dhdGetButton(0);
		if (button0_state_)
		{
			//获得omega的新位姿
			dhdGetPosition(&current_position[0], &current_position[1], &current_position[2]);
			dhdGetOrientationFrame(Omega_matrix);

			//计算位移
			TcpPose[0] = TcpPoseInit[0] + (current_position[0] - current_position_Init[0]) * PositonScale;
			TcpPose[1] = TcpPoseInit[1] + (current_position[1] - current_position_Init[1]) * PositonScale;
			TcpPose[2] = TcpPoseInit[2] + (current_position[2] - current_position_Init[2]) * PositonScale;
			// outfile3 << TcpPose[0] << " " << TcpPose[1] << " " << TcpPose[2] << std::endl;
			//计算姿态
			//首先计算deltaR，deltaR=omega新位姿*omegaInit位姿的转置
			Omega_matrixEigen << Omega_matrix[0][0], Omega_matrix[0][1], Omega_matrix[0][2],
				Omega_matrix[1][0], Omega_matrix[1][1], Omega_matrix[1][2],
				Omega_matrix[2][0], Omega_matrix[2][1], Omega_matrix[2][2];

			DeltaR = Omega_matrixEigen * Omega_matrixInitEigen.transpose();
			//新的tcp姿态=deltaR*tcpinit姿态,进而得到TcpPose[3],TcpPose[4],TcpPose[5]
			Tcp_Rotation = DeltaR * Tcp_Rotation_Init;
			R = Tcp_Rotation_Init.inverse() * Tcp_Rotation;
			K.fromRotationMatrix(R);
			Khat = K.axis();
			alpha = K.angle();
			R = Eigen::AngleAxisd(alpha * PoseScale, Khat);
			Tcp_Rotation = Tcp_Rotation_Init * R;
			K.fromRotationMatrix(Tcp_Rotation);
			Khat = K.axis();
			alpha = K.angle();
			TcpPose[3] = Khat[0] * alpha;
			TcpPose[4] = Khat[1] * alpha;
			TcpPose[5] = Khat[2] * alpha;

			//将TcpPose发给机器人，机器人开始运动
			// dhdSetForce(fx, fy, fz, -1);
		}
		else
		{
			// rtde_control.zeroFtSensor();
			//  rtde_control.zeroFtSensor();
			//获得Init数据，机器人TCP和主手位姿
			dhdGetPosition(&current_position_Init[0], &current_position_Init[1], &current_position_Init[2]);
			// outfile4 << current_position_Init[0] << " " << current_position_Init[1] << " " << current_position_Init[2] << std::endl;

			// TcpPoseInit = rtde_receive.getActualTCPPose();
			TcpPoseInit = TcpPose;
			dhdGetOrientationFrame(Omega_matrixInit);
			Omega_matrixInitEigen << Omega_matrixInit[0][0], Omega_matrixInit[0][1], Omega_matrixInit[0][2],
				Omega_matrixInit[1][0], Omega_matrixInit[1][1], Omega_matrixInit[1][2],
				Omega_matrixInit[2][0], Omega_matrixInit[2][1], Omega_matrixInit[2][2];

			//计算机器人Init姿态
			axisInit = {TcpPoseInit[3], TcpPoseInit[4], TcpPoseInit[5]};
			axisInit.normalize();
			alphaInit = TcpPoseInit[3] / axisInit[0];
			Tcp_Rotation_Init = Eigen::AngleAxisd(alphaInit, axisInit);
		}

		TcpForce = rtde_receive.getActualTCPForce();

		double angle = sqrt(KDL::pow(TcpPose[3], 2) + pow(TcpPose[4], 2) + pow(TcpPose[5], 2));
		KDL::Frame TcpPoseKDL{KDL::Rotation::Rot(KDL::Vector{TcpPose[3], TcpPose[4], TcpPose[5]}, angle), KDL::Vector{TcpPose[0], TcpPose[1], TcpPose[2]}};

		ur_admittance.set_force(TcpForce[0], TcpForce[1], TcpForce[2]);
		ur_admittance.set_torque(TcpForce[3], TcpForce[4], TcpForce[5]);
		ur_admittance.calculate(frame_offset, admittance_vel);

		KDL::Frame KDL_TcpPoseCommand = frame_offset * TcpPoseKDL;

		std_tcp_pos_command[0] = KDL_TcpPoseCommand.p(0);
		std_tcp_pos_command[1] = KDL_TcpPoseCommand.p(1);
		std_tcp_pos_command[2] = KDL_TcpPoseCommand.p(2);
		KDL::Vector KhatKdl = KDL_TcpPoseCommand.M.GetRot();
		std_tcp_pos_command[3] = KhatKdl[0];
		std_tcp_pos_command[4] = KhatKdl[1];
		std_tcp_pos_command[5] = KhatKdl[2];

		rtde_control.servoL(std_tcp_pos_command, 0.5, 0.5, 0.002, 0.05, 300);
		std::this_thread::sleep_for(std::chrono::duration<double>(0.002));
		fx = f1.filter(TcpForce[0]) * scaling;
		fy = f2.filter(TcpForce[1]) * scaling;
		fz = f3.filter(TcpForce[2]) * scaling;
		if (abs(fx) > 0.6)
			fx = sign(fx)*0.6;
		if (abs(fy) > 0.6)
			fy =sign(fy)*0.6;
		if (abs(fz) > 0.6)
			fz = sign(fz)*0.6;

		// if (KDL::Vector{fx, fy, fz}.Norm() > 5)
		// {
		// 	fx = sign(fx) * 5;
		// 	fy = sign(fy) * 5;
		// 	fz = sign(fz) * 5;
		// }

		//dhdSetForce(fx, fy, fz, -1);
	}
	// process finished.
	haptic_dev.keep_alive_ = false;
	spinner.stop();
#endif

#if 0 

	KDL::Frame start;
	start.p[0] = TcpPose[0];
	start.p[1] = TcpPose[1];
	start.p[2] = TcpPose[2];
	double angle = sqrt(KDL::pow(TcpPose[3], 2) + pow(TcpPose[4], 2) + pow(TcpPose[5], 2));
	start.M = KDL::Rotation::Rot(KDL::Vector{TcpPose[3], TcpPose[4], TcpPose[5]}, angle);
	while (ros::ok())
	{

		ur_admittance.calculate(frame_offset, admittance_vel);

		KDL::Frame frame_target;

		frame_target = frame_offset * start;

		std::vector<double> tcp_force_torque = rtde_receive.getActualTCPForce();

		ur_admittance.set_force(tcp_force_torque[0], tcp_force_torque[1], tcp_force_torque[2]);
		ur_admittance.set_torque(tcp_force_torque[3], tcp_force_torque[4], tcp_force_torque[5]);

		TcpPose[0] = frame_target.p(0);
		TcpPose[1] = frame_target.p(1);
		TcpPose[2] = frame_target.p(2);
		KDL::Vector KhatKdl = frame_target.M.GetRot();
		TcpPose[3] = KhatKdl[0];
		TcpPose[4] = KhatKdl[1];
		TcpPose[5] = KhatKdl[2];

		rtde_control.servoL(TcpPose, 0.5, 0.5, 0.002, 0.05, 300);
		std::this_thread::sleep_for(std::chrono::duration<double>(0.002));
	}

#endif

	return 0;
}
