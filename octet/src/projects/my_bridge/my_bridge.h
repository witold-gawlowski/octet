namespace octet {
  class my_bridge : public app {
    ref<visual_scene> app_scene;
    scene_node *camera_node;
    mouse_look mouse_look_helper;
    collada_builder loader;
  public:
    my_bridge(int argc, char **argv) : app(argc, argv) {
    }

    ~my_bridge() {
    }
    void app_init() {
      app_scene =  new visual_scene();

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

      //setting up camera and mouse helper
      mouse_look_helper.init(this, 50.0f / 360.0f, false);
      app_scene->create_default_camera_and_lights();
      camera_node = app_scene->get_camera_instance(0)->get_node();
      camera_node->loadIdentity();

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
        mat.translate(pos[0], pos[1], pos[2]);
        mesh_instance *col = app_scene->add_shape(mat, new mesh_box(vec3(scl[0], scl[1], scl[2])), red, false);
        
        //std::cout << pos[0] << " " << rot[1] << " " << scl[2] << std::endl;
      }

    }



    void draw_world(int x, int y, int w, int h) {
      if (this->is_key_down('W')) {
        camera_node->translate(vec3(0, 0, -1));
      }
      else if (this->is_key_down('S')) {
        camera_node->translate(vec3(0, 0, 1));
      }
      else if (this->is_key_down('A')) {
        camera_node->translate(vec3(-1, 0, 0));
      }
      else if (this->is_key_down('D')) {
        camera_node->translate(vec3(1, 0, 0));
      }
      else if (this->is_key_down('Q')) {
        camera_node->translate(vec3(0, -1, 0));
      }
      else if (this->is_key_down('E')) {
        camera_node->translate(vec3(0, 1, 0));
      }
      else if (this->is_key_down(' ')) {
        //dupm position into the console
        printf("%f %f %f\n", camera_node->get_position().x(), camera_node->get_position().y(), camera_node->get_position().z());
      }
      mat4t &camera_to_world = camera_node->access_nodeToParent();
      mouse_look_helper.update(camera_to_world);
      int vx = 0, vy = 0;
      get_viewport_size(vx, vy);
      app_scene->begin_render(vx, vy);
      app_scene->update(1.0f/30);
      app_scene->render((float)vx / vy);
    }
  };
}
