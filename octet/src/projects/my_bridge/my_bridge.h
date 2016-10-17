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
      mouse_look_helper.init(this, 50.0f / 360.0f, false);
      

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
        //bellow transformation needs to be done for any data that come from blender to opengl. 
        node->rotate(-90, vec3(1, 0, 0));
        app_scene->add_child(node);
        app_scene->add_mesh_instance(new mesh_instance(node, canyon, mat));
      }
      app_scene->create_default_camera_and_lights();
      camera_node = app_scene->get_camera_instance(0)->get_node();
    
      camera_node->loadIdentity();
      //camera_node->rotate(90, vec3(1, 0, 0));
      //camera_node->rotate(90, vec3(1, 0, 0));
      //camera_node->rotate(30, vec3(1, 0, 0));
      //camera_node->rotate(30, vec3(0, 1, 0));
      //camera_node->translate(vec3(-500, 1000, 500));
      //camera_node->translate(vec3(0, 0, 500));
      //camera_node->rotate(10, vec3(1, 0, ));
      //camera_node->rotate(90, vec3(0, 1, 0));
      //camera_node->translate(vec3(0, 0, 100));
      /*the internal units in octet are centimeters, right?
      I have imported a model in collada, and everything looks like they are.*/



      vec3 camera_position = app_scene->get_camera_instance(0)->get_node()->get_position();
      printf("%f %f %f", camera_position.x(), camera_position.y(), camera_position.z());
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
        camera_node->translate(vec3(0, -1, 0));
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
