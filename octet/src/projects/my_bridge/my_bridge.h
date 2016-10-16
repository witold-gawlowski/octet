namespace octet {
  class my_bridge : public app {
    ref<visual_scene> app_scene;

  collada_builder loader;
  public:
    my_bridge(int argc, char **argv) : app(argc, argv) {
    }

    ~my_bridge() {
    }
    void app_init() {
      app_scene =  new visual_scene();
      
      

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
        node->translate(vec3(0, 0, 0));
        app_scene->add_child(node);
        app_scene->add_mesh_instance(new mesh_instance(node, canyon, mat));
      }
      app_scene->create_default_camera_and_lights();
      scene_node *camera_node = app_scene->get_camera_instance(0)->get_node();
      camera_node->loadIdentity();
    



      vec3 camera_position = app_scene->get_camera_instance(0)->get_node()->get_position();
      printf("%f %f %f", camera_position.x(), camera_position.y(), camera_position.z());
    }

    void draw_world(int x, int y, int w, int h) {
      int vx = 0, vy = 0;
      get_viewport_size(vx, vy);
      app_scene->begin_render(vx, vy);
      app_scene->update(1.0f/30);
      app_scene->render((float)vx / vy);
    }
  };
}