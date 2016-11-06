namespace octet {

  class my_bridge : public app {
    static vec3 pos (btRigidBody* rb) {
      btVector3 position = rb->getCenterOfMassPosition ();
      return vec3 (position.x (), position.y (), position.z ());
    }
    class segment {
      //delarations at the bottom of this file
      const static float bridge_width;
      const static float bridge_height;
      const static float segment_length;
    public:
      struct linker {
        btRigidBody *ul, *ur, *dl, *dr;
        linker () { ul = ur = dl = dr = nullptr; }
        linker (btRigidBody *ul, btRigidBody *ur, btRigidBody *dl, btRigidBody *dr) : ul (ul), ur(ur), dl(dl), dr(dr){}
        
      };

    private:
      linker forward, backward;
      vec3 dir;
      ref<visual_scene> app_scene;

    public:
      explicit segment (linker l, vec3 dir, visual_scene* scene) : backward(l), dir(dir), app_scene(scene) {
        material *red = new material (vec4 (1, 0, 0, 1));
        mesh_instance *plank = app_scene->add_shape (mat4t().translate(pos(l.dl)+vec3(bridge_width/2.0f, 0, -segment_length*(0.25+1/16.))),
          new mesh_box (vec3 (bridge_width/2.0f, 0.01f, segment_length/8)), red,false);
      }

    };

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
    

   void make_join(vec3 pos){
      
    }
    
    

    void app_init() {
      app_scene = new visual_scene();

      //importing canyon:
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

      //importing poles
      
      if ( !loader.load_xml ("assets/projects/my_bridge/meshes/poles.dae") ) {
        return;
      }
      dict.reset ();
      meshes.reset ();
      loader.get_resources (dict);
      dict.find_all (meshes, atom_mesh);
      if ( meshes.size () ) {
        material *mat = new material (vec4 (0.25f, 0.25f, 0.5f, 1));
        mesh *pole = meshes[0]->get_mesh ();
        scene_node *node = new scene_node ();
        //bellow transformation needs to be done for any data that comes from blender to opengl. 
        node->rotate (-90, vec3 (1, 0, 0));
        app_scene->add_child (node);
        app_scene->add_mesh_instance (new mesh_instance (node, pole, mat));
      }



      //setting up camera, mouse helper and fps helper
      mouse_look_helper.init(this, 35.0f / 360.0f, false);
      fps_helper.init(this);
      app_scene->create_default_camera_and_lights();
      the_camera = app_scene->get_camera_instance(0);
      the_camera->get_node()->loadIdentity();
      the_camera->set_far_plane(10000);
      the_camera->set_near_plane(0.3f);
      



      //preparing player
      float player_height = 0.0f;
      float player_mass = 0.5f;
      mat4t player_mat;
      player_mat.loadIdentity();
      player_mat.translate(0, 5, 16);

      mesh_instance *mi = app_scene->add_shape(
        player_mat,
        new mesh_sphere(vec3(0), 0),
        new material(vec4(0, 0, 1, 1)),
        true, player_mass,
        new btCapsuleShape(0.05f, player_height)
      );
      player_node = mi->get_node();

      //read and add world colliders
      material *collider_mat = new material(vec4(0.15, 0.15, 0.25, 1));
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
        mesh_instance *col = app_scene->add_shape(mat, new mesh_box(vec3(scl[0], scl[1], scl[2])), collider_mat, false);
        //std::cout << pos[0] << " " << rot[1] << " " << scl[2] << std::endl;
      }

      //generate bridge
      //mesh_instance *col = app_scene->add_shape (mat4t().translate(vec3(-1, -5, -7)), new mesh_box (vec3 (2)), red,false);
      //player-side dockers
      vec3 tld_pos (0.44f, 0.31, 12.85), trd_pos (1.31f, 0.31, 12.85), drd_pos (1.31f, -0.31, 12.85), dld_pos (0.44f, -0.31, 12.85);
      material *red = new material (vec4 (1, 0, 0, 1));
      mesh_instance *tld = app_scene->add_shape (mat4t ().translate (tld_pos), new mesh_box (vec3 (0.07)), red, false);
      mesh_instance *trd = app_scene->add_shape (mat4t ().translate (trd_pos), new mesh_box (vec3 (0.07)), red, false);
      mesh_instance *drd = app_scene->add_shape (mat4t ().translate (drd_pos), new mesh_box (vec3 (0.07)), red, false);
      mesh_instance *dld = app_scene->add_shape (mat4t ().translate (dld_pos), new mesh_box (vec3 (0.07)), red, false);
      segment *seg = new segment(segment::linker(tld->get_node()->get_rigid_body(), 
        trd->get_node ()->get_rigid_body (),
          dld->get_node ()->get_rigid_body (),
            drd->get_node ()->get_rigid_body ()),
              vec3(0, 0, -1), app_scene);
      

      //now read, add and bind bridge elements3 
    //  btVector3 axisA(1.0f, 0.0f, 0.0f);
    //  btVector3 axisB(1.0f, 0.0f, 0.0f);
    //  btVector3 pivotA(0, -0.68f, 0.f);
    //  btVector3 pivotB(0, 0.68f, 0.f);
    //  btHingeConstraint *bridge_hinge;
    //  mesh_instance *previous_col = NULL;
    //  ifs >> n;
    //  for (int i = 0; i < n; i++) {
    //    ifs >> pos[0] >> pos[1] >> pos[2];
    //    ifs >> rot[0] >> rot[1] >> rot[2] >> rot[3];
    //    ifs >> scl[0] >> scl[1] >> scl[2];
    //    mat4t mat(rot);
    //    mat.rotateX(90);
    //    mat.translate(pos[0], pos[1], pos[2]);
    //    mesh_instance *col = app_scene->add_shape(mat, new mesh_box(vec3(scl[0], scl[1], scl[2])), red, (i==0||i==n-1) ? false : true, 5.0f);
    //    btRigidBody *rbA; 
    //    btRigidBody *rbB = col->get_node()->get_rigid_body();
    //    rbB->setDamping(0.1f, 0.1f);
    //    rbB->applyDamping(0.15f);
    //    //printf("%d bool: %d\n", previous_col, previous_col != NULL);
    //    if (previous_col) {
    //      rbA = previous_col->get_node()->get_rigid_body();
    //      bridge_hinge = new btHingeConstraint(*rbA, *rbB, pivotA, pivotB, axisA, axisB);
    //      bridge_hinge->setLimit(-SIMD_HALF_PI * 0.5f, SIMD_HALF_PI * 0.5f);
    //      //what is CRP and ERP here: http://bulletphysics.org/mediawiki-1.5.8/index.php/Definitions
    //      bridge_hinge->setParam(BT_CONSTRAINT_STOP_ERP, 0.8);
    //      bridge_hinge->setParam(BT_CONSTRAINT_STOP_CFM, 0.1);  
    //      bridge_hinge->setParam(BT_CONSTRAINT_CFM, 0.05);
    //      app_scene->add_my_hinge(bridge_hinge);
    //    }
    //    previous_col = col;
    //  }
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
  //its the aproximate distance between the poles
  const float my_bridge::segment::bridge_width = 1.31f -0.44f;
  const float my_bridge::segment::bridge_height = bridge_width * 1.618;
  const float my_bridge::segment::segment_length = bridge_width;
}
