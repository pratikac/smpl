#include <iostream>
#include <ctime>

#include <ros/ros.h>
#include <message_filters/subscriber.h>
#include <geometry_msgs/Point32.h>
#include <geometry_msgs/PointStamped.h>
#include <geometry_msgs/PoseStamped.h>
#include <geometry_msgs/Pose.h>
#include <sensor_msgs/PointCloud.h>
#include <nav_msgs/OccupancyGrid.h>
#include <nav_msgs/Path.h>

#include <tf/tf.h>

using namespace std;


class SimpleGoal
{
  public:

    ros::Subscriber rviz_goal_listener;
    ros::Publisher goal_pub;
    SimpleGoal()
    {
      ros::NodeHandle n;
      rviz_goal_listener = n.subscribe("move_base_simple/goal", 2, &SimpleGoal::on_rviz_goal_listener, this);
      goal_pub = n.advertise<geometry_msgs::PoseStamped>("pnc_nextpose", 2);
    }
    void on_rviz_goal_listener(const geometry_msgs::PoseStamped::ConstPtr p)
    {
      geometry_msgs::PoseStamped pose;
      pose = *p;
      goal_pub.publish(pose);
      ROS_INFO("published pose");
      cout<<pose<<endl;
    }
};




int main(int argc, char **argv)
{
  ros::init(argc, argv, "simple_goal");

  SimpleGoal sg;

  ros::spin();

  return 0;
}
