#include <gtest/gtest.h>

#include <wmrde/wmrstate.h>
#include <wmrde/demo/zoemodel.h>
#include <wmrde/util/string_format.h>

//uncomment one to enable/disable print output
#define DEBUG_INFO(x) do { std::cout << x << std::endl; } while(0)
//#define DEBUG_INFO(x) do { } while(0)

using namespace wmrde;

TEST(TestSuite, convert_state_to_HT)
{
  //init WmrModel
  WmrModel mdl;
  makeZoeModel(mdl);

  DEBUG_INFO(string_format("mdl. numFrames = %d, numWheelFrames = %d, numActuatedFrames = %d",
      mdl.numFrames(), mdl.numWheelFrames(), mdl.numActuatedFrames()));

  //init WmrState
  WmrState state;
  state.position = (Vec3() << 0, 5, 1).finished();
  Mat3 R = eulerToRot(0, 0, degToRad(15.0));
  state.orientation = Quaternion(R);
  state.joint_disp.resize(mdl.numJoints());
  state.joint_disp.setZero();

  state.joint_disp[mdl.frameNameToIndex("front_steer")-1] = degToRad(-15.0);
  state.joint_disp[mdl.frameNameToIndex("front_left_wheel")-1] = degToRad(20.0);

  DEBUG_INFO("state.orientation.coeffs() = \n" << state.orientation.coeffs());
  DEBUG_INFO("state.orientation.matrix() = \n" << state.orientation.matrix());
  DEBUG_INFO("Rotation matrix input to orientation ctor = \n" << R);

  //compute HT_parent
  std::vector<HTransform> HT_parent;
  stateToHTparent(mdl, state, HT_parent);

  for (size_t i = 0; i < HT_parent.size(); i++)
  {
    DEBUG_INFO("HT_parent[" << mdl.getFrame(i).name << "] = \n" << HT_parent[i]);
  }

  //compute HT_world
  std::vector<HTransform> HT_world;
  HTparentToHTworld(mdl, HT_parent, HT_world);
  for (size_t i = 0; i < HT_world.size(); i++)
  {
    DEBUG_INFO("HT_world[" << mdl.getFrame(i).name << "] = \n" << HT_world[i]);
  }

  //TODO, check this
  Vec3 expected = (Vec3() << 0.922459, 4.42717, 0.881).finished();
  EXPECT_TRUE(HT_world[5].t.isApprox(expected, 1e-4));
}

// Run all the tests that were declared with TEST()
int main(int argc, char **argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
