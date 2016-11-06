namespace octet {
  class my_bridge : public app {
    ref<visual_scene> app_scene;
    ref<camera_instance> the_camera;
    mouse_look mouse_look_helper;
    ref<scene_node> player_node;
    helper_fps_controller fps_helper;
    collada_builder loader; 

  public:
    my_bridge(int argc, char **argv) : app(argc, argv) {
    }

    ~my_bridge() {
    }
    void app_init() {
      app_scene = new visual_scene();

      //importing collada canyon:
      resource_dict dict;
      if (!loader.load_xml("assets/projects/my_bridge/meshes/Canyon.dae")) {
        return;
      }
      loader.get_resources(dict);
      dynarray<resource*> meshes;
      dict.find_all(meshes, atom_mesh);
      if (meshes.size()) {
        material *mat = new material(vec4(0.25f, 0.25f, 0.5f, 1));
        mesh *canyon = meshes[0]->get_mesh();
        scene_node *node = new scene_node();
        //bellow transformation needs to be done for any data that comes from blender to opengl. 
        node->rotate(-90, vec3(1, 0, 0));
        app_scene->add_child(node);
        app_scene->add_mesh_instance(new mesh_instance(node, canyon, mat));
      }

      //setting up camera, mouse helper and fps helper
      mouse_look_helper.init(this, 50.0f / 360.0f, false);
      fps_helper.init(this);
      app_scene->create_default_camera_and_lights();
      the_camera = app_scene->get_camera_instance(0);
      the_camera->get_node()->loadIdentity();
      the_camera->set_far_plane(10000);
      the_camera->set_near_plane(0.4f);


      //preparing player
      float player_height = 1.83f;
      float player_radius = 0.35f;
      float player_mass = 10.0f;
      mat4t player_mat;
      player_mat.loadIdentity();
      player_mat.translate(0, 5, 16);

      mesh_instance *mi = app_scene->add_shape(
        player_mat,
        new mesh_sphere(vec3(0), player_radius),
        new material(vec4(0, 0, 1, 1)),
        true, player_mass,
        new btCapsuleShape(0.95f, player_height)
      );
      player_node = mi->get_node();

      //read and add world colliders
      material *red = new material(vec4(1, 0, 0, 1));
      std::ifstream ifs("box_list.csv");
      int n;
      //rotation via quaternion 
      vec4 pos, rot, scl;
      ifs >> n;
      for (int i = 0; i < n; i++) {
        ifs >> pos[0] >> pos[1] >> pos[2];
        ifs >> rot[0] >> rot[1] >> rot[2] >> rot[3];
        ifs >> scl[0] >> scl[1] >> scl[2];
        mat4t mat(rot);
        mat.rotateX(90);
        mat.translate(pos[0], pos[1], pos[2]);
        mesh_instance *col = app_scene->add_shape(mat, new mesh_box(vec3(scl[0], scl[1], scl[2])), red, false);
        //std::cout << pos[0] << " " << rot[1] << " " << scl[2] << std::endl;
      }

      //now read, add and bind bridge elements3
      btVector3 axisA(1.0f, 0.0f, 0.0f);
      btVector3 axisB(1.0f, 0.0f, 0.0f);
      btVector3 pivotA(0, -0.68f, 0.f);
      btVector3 pivotB(0, 0.68f, 0.f);
      btHingeConstraint *bridge_hinge;
      mesh_instance *previous_col = NULL;
      ifs >> n;
      for (int i = 0; i < n; i++) {
        ifs >> pos[0] >> pos[1] >> pos[2];
        ifs >> rot[0] >> rot[1] >> rot[2] >> rot[3];
        ifs >> scl[0] >> scl[1] >> scl[2];
        mat4t mat(rot);
        mat.rotateX(90);
        mat.translate(pos[0], pos[1], pos[2]);
        mesh_instance *col = app_scene->add_shape(mat, new mesh_box(vec3(scl[0], scl[1], scl[2])), red, (i==0||i==n-1) ? false : true, 5.0f);
        btRigidBody *rbA; 
        btRigidBody *rbB = col->get_node()->get_rigid_body();
        rbB->setDamping(0.1f, 0.1f);
        rbB->applyDamping(0.15f);
        //printf("%d bool: %d\n", previous_col, previous_col != NULL);
        if (previous_col) {
          rbA = previous_col->get_node()->get_rigid_body();
          bridge_hinge = new btHingeConstraint(*rbA, *rbB, pivotA, pivotB, axisA, axisB);
          bridge_hinge->setLimit(-SIMD_HALF_PI * 0.5f, SIMD_HALF_PI * 0.5f);
          //what is CRP and ERP here: http://bulletphysics.org/mediawiki-1.5.8/index.php/Definitions
          bridge_hinge->setParam(BT_CONSTRAINT_STOP_ERP, 0.8);
          bridge_hinge->setParam(BT_CONSTRAINT_STOP_CFM, 0.1);  
          bridge_hinge->setParam(BT_CONSTRAINT_CFM, 0.05);
          app_scene->add_my_hinge(bridge_hinge);
        }
        previous_col = col;
      }
    }


    void draw_world(int x, int y, int w, int h) {
      //moving around the ghost
      /*
      if (this->is_key_down('W')) {
        the_camera->get_node()->translate(vec3(0, 0, -1));
      }
      else if (this->is_key_down('S')) {
        the_camera->get_node()->translate(vec3(0, 0, 1));
      }
      else if (this->is_key_down('A')) {
        the_camera->get_node()->translate(vec3(-1, 0, 0));
      }
      else if (this->is_key_down('D')) {
        the_camera->get_node()->translate(vec3(1, 0, 0));
      }
      else if (this->is_key_down('Q')) {
        the_camera->get_node()->translate(vec3(0, -1, 0));
      }
      else if (this->is_key_down('E')) {
        the_camera->get_node()->translate(vec3(0, 1, 0));
      }
      else if (this->is_key_down(' ')) {
        //dupm position into the console
        printf("%f %f %f\n", the_camera->get_node()->get_position().x(), the_camera->get_node()->get_position().y(), the_camera->get_node()->get_position().z());
      }
      */
      mat4t &camera_to_world = the_camera->get_node()->access_nodeToParent();
      mouse_look_helper.update(camera_to_world);
      fps_helper.update(player_node, the_camera->get_node());
      int vx = 0, vy = 0;
      get_viewport_size(vx, vy);
      app_scene->begin_render(vx, vy);
      app_scene->update(1.0f/30);
      app_scene->render((float)vx / vy);
    }
  };
}
