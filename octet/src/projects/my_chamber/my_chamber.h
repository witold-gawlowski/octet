#include <memory>
namespace octet {
  class my_chamber : public app {
    ref<visual_scene> app_scene;
    static const int grid_size = 100;
    static const int fountain_x = 50;
    static const int fountain_y = 30;
    class mesh_fluid : public mesh {
      struct my_vertex {
        vec3p pos;
        vec3p color;
      };

      //boxes can only be square
      dynarray<my_vertex> vertices;
      std::vector<float> prev_density;
      std::vector<float> prev_vx;
      std::vector<float> prev_vy;
      std::vector<float> density;
      std::vector<float> vx;
      std::vector<float> vy;
      std::vector<int> mask;
      std::vector<vec3> boxes; //(pos_x, pos_y, size)
      ivec3 dim;
      float cx, cy, sx, sy;

    public:
      mesh_fluid (aabb_in bb, ivec3_in dim) : mesh (), dim (dim) {
        mesh::set_aabb (bb);
        sx = bb.get_half_extent ().x ()*(2.0f / dim.x ());
        sy = bb.get_half_extent ().y ()*(2.0f / dim.y ());
        cx = bb.get_center ().x () - bb.get_half_extent ().x ();
        cy = bb.get_center ().y () - bb.get_half_extent ().y ();

        density.resize ((dim.x () + 1)*(dim.y () + 1));
        vx.resize ((dim.x () + 1)*(dim.y () + 1));
        vy.resize ((dim.x () + 1)*(dim.y () + 1));
        mask.resize ((dim.x () + 1)*(dim.y () + 1));

        prev_density.resize ((dim.x () + 1)*(dim.y () + 1));
        prev_vx.resize ((dim.x () + 1)*(dim.y () + 1));
        prev_vy.resize ((dim.x () + 1)*(dim.y () + 1));

        //density[50 +(dim.x()+1) * 50] = 1;
        //prev_vx[50 +(dim.x()+1) * 50] = 1;

        dynarray<uint32_t> indices;
        int stride = dim.x () + 1;
        for ( int i = 0; i < dim.x (); ++i ) {
          for ( int j = 0; j < dim.y (); ++j ) {
            indices.push_back ((i + 1) + (j + 0)*stride);
            indices.push_back ((i + 0) + (j + 1)*stride);
            indices.push_back ((i + 1) + (j + 1)*stride);
            indices.push_back ((i + 1) + (j + 0)*stride);
            indices.push_back ((i + 0) + (j + 0)*stride);
            indices.push_back ((i + 0) + (j + 1)*stride);
          }
        }
        set_indices (indices);
        clear_attributes ();
        add_attribute (attribute_pos, 3, GL_FLOAT, 0);
        add_attribute (attribute_color, 3, GL_FLOAT, 12);
      }

      bool movable(int N, const vec3 &box, int dx, int dy){
        auto IX = [=] (int i, int j) { return i + (N + 2)*j; };
        assert (!dx || !dy);
        if ( dx>0 ) {
          for ( int i = box.y (); i<box.y () + box.z (); ++i ) {
            if ( mask[IX (box.x ()+box.z(), i)] ) {
              return false;
            }
          }
          return true;
        }else if(dx<0){
          for ( int i = box.y (); i<box.y () + box.z (); ++i ) {
            if ( mask[IX (box.x ()-1, i)] ) {
              return false;
            }
          }
          return true;
        } else if ( dy>0 ) {
          for ( int i = box.x (); i<box.x () + box.z (); ++i ) {
            if ( mask[IX (i, box.y()+box.z())] ) {
              return false;
            }
          }
          return true;
        } else if ( dy<0 ) {
          for ( int i = box.x (); i<box.x () + box.z (); ++i ) {
            if ( mask[IX (i, box.y()-1)] ) {
              return false;
            }
          }
          return true;
        }
      }

      void move_box(int N, int index, int dx, int dy){
        //todo: repair infinite amplification on the edges parallell to the movement
        vec3 &box = boxes[index];
        assert (!dx || !dy);
        auto IX = [=] (int i, int j) { return i + (N + 2)*j; };
        if ( !movable (N, box, dx, dy) )
          return;
        for(int i=box.x(); i<box.x()+box.z(); ++i){
          for(int j=box.y(); j<box.y()+box.z(); ++j){
            if( mask[IX (i+dx, j+dy)]!=1){
              density[IX (i + 2*dx, j + 2*dy)] += density[IX (i + dx, j + dy)]/1.5f;
            }
            density[IX (i, j)] = 0;
            mask[IX (i, j)] = 0;
          }
        }
        box[0] += dx;
        box[1] += dy;
        for ( int i = box.x (); i < box.x () + box.z (); ++i ) {
          for ( int j = box.y (); j < box.y () + box.z (); ++j ){
            mask[IX (i, j)] = 1;
          }
        }
        for ( int i = box.x (); i<box.x () + box.z (); ++i ) {
          for ( int j = box.y (); j<box.y () + box.z (); ++j ) {            
            if ( mask[IX (i + dx, j + dy)] != 1 ) {
              if ( dx ) {
                
                vx[IX (i + dx, j + dy)] += 3*dx;
              } else {
                vy[IX (i + dx, j + dy)] += 3 * dy;
              }
            }else if( mask[IX (i - dx, j - dy)] != 1 ){
              density[IX (i - dx, j - dy)] = 0;
              if ( dx ) {
                std::cout << "ads" << std::endl;
                vx[IX (i - dx, j - dy)] += 2 * dx;
                //vx[IX (i - 2*dx, j - 2*dy)] += 1.5 * dx;
              } else {
                vy[IX (i - dx, j - dy)] += 2 * dy;
                //vy[IX (i - 2*dx, j - 2*dy)] += 1.5 * dy;
              }
            }
          }
        }
      }

      void add_box (int N, int x, int y, int size) {
        auto IX = [=] (int i, int j) { return i + (N + 2)*j; };
        for ( int i = x; i<x + size; ++i )
          for ( int j = y; j<y + size; ++j ) {
            mask[IX (i, j)] = 1;
          }
        boxes.push_back (vec3 (x, y, size));
      }


      void set_my_boundary2 (int N, int b, float * x, vec2 pos, int size) {
        auto IX = [=] (int i, int j) { return i + (N + 2)*j; };
        for ( int i = 0; i<size; i++ ) {
          x[IX (pos.x () - 1, pos.y () + i)] = b == 1 ? -x[IX (pos.x () - 2, pos.y () + i)] : x[IX (pos.x () - 2, pos.y () + i)];
          x[IX (pos.x () + size, pos.y () + i)] = b == 1 ? -x[IX (pos.x () + size + 1, pos.y () + i)] : x[IX (pos.x () + size + 1, pos.y () + i)];
          x[IX (pos.x () + i, pos.y () - 1)] = b == 2 ? -x[IX (pos.x () + i, pos.y () - 2)] : x[IX (pos.x () + i, pos.y () - 2)];
          x[IX (pos.x () + i, pos.y () + size)] = b == 2 ? -x[IX (pos.x () + i, pos.y () + size + 1)] : x[IX (pos.x () + i, pos.y () + size + 1)];
        }
        x[IX (pos.x () - 1, pos.y () - 1)] = 0.5f*(x[IX (pos.x (), pos.y () - 1)] + x[IX (pos.x () - 1, pos.y ())]);
        x[IX (pos.x () + size, pos.y () - 1)] = 0.5f*(x[IX (pos.x () + size, pos.y () - 2)] + x[IX (pos.x () + size + 1, pos.y () - 1)]);
        x[IX (pos.x () + size, pos.y () + size)] = 0.5f*(x[IX (pos.x () + size + 1, pos.y () + size)] + x[IX (pos.x () + size, pos.y () + size + 1)]);
        x[IX (pos.x () - 1, pos.y () + size)] = 0.5f*(x[IX (pos.x () - 2, pos.y () + size)] + x[IX (pos.x () - 1, pos.y () + size + 1)]);
      }

      void set_my_boundary (int N, int b, float * x, vec2 pos, int size) {
        auto IX = [=] (int i, int j) { return i + (N + 2)*j; };
        for ( int i = 1; i<size-1; i++ ) {
          x[IX (pos.x (), pos.y () + i)] = b == 1 ? -x[IX (pos.x () - 1, pos.y () + i)] : x[IX (pos.x () - 1, pos.y () + i)];
          x[IX (pos.x () + size-1, pos.y () + i)] = b == 1 ? -x[IX (pos.x () + size, pos.y () + i)] : x[IX (pos.x () + size , pos.y () + i)];
          x[IX (pos.x () + i, pos.y ())] = b == 2 ? -x[IX (pos.x () + i, pos.y () - 1)] : x[IX (pos.x () + i, pos.y () - 1)];
          x[IX (pos.x () + i, pos.y () + size-1)] = b == 2 ? -x[IX (pos.x () + i, pos.y () + size )] : x[IX (pos.x () + i, pos.y () + size)];
        }
        x[IX (pos.x (), pos.y ())] = 0.5f*(x[IX (pos.x (), pos.y ()+1)] + x[IX (pos.x () + 1, pos.y ())]);
        x[IX (pos.x () + size-1, pos.y ())] = 0.5f*(x[IX (pos.x () + size-2, pos.y ())] + x[IX (pos.x () + size -1, pos.y () + 1)]);
        x[IX (pos.x () + size-1, pos.y () + size-1)] = 0.5f*(x[IX (pos.x () + size -2, pos.y () + size-1)] + x[IX (pos.x () + size-1, pos.y () + size - 2)]);
        x[IX (pos.x (), pos.y () + size - 1)] = 0.5f*(x[IX (pos.x (), pos.y () + size-2)] + x[IX (pos.x () + 1, pos.y () + size -1)]);
      }

      void set_chamber_walls (int N, int b, float * x) {
        auto IX = [=] (int i, int j) { return i + (N + 2)*j; };

        //W: i need to keep distance from boxes so that moving them does not consume fluid
        for ( int i = 1; i <= N; i++ ) {
          x[IX (0, i)] = b == 1 ? -x[IX (1, i)] : x[IX (1, i)];
          x[IX (N + 1, i)] = b == 1 ? -x[IX (N, i)] : x[IX (N, i)];
          x[IX (i, 0)] = b == 2 ? -x[IX (i, 1)] : x[IX (i, 1)];
          x[IX (i, N + 1)] = b == 2 ? -x[IX (i, N)] : x[IX (i, N)];
        }

        x[IX (0, 0)] = 0.5f*(x[IX (1, 0)] + x[IX (0, 1)]);
        x[IX (0, N + 1)] = 0.5f*(x[IX (1, N + 1)] + x[IX (0, N)]);
        x[IX (N + 1, 0)] = 0.5f*(x[IX (N, 0)] + x[IX (N + 1, 1)]);
        x[IX (N + 1, N + 1)] = 0.5f*(x[IX (N, N + 1)] + x[IX (N + 1, N)]);
      }
      void set_boundary (int N, int b, float * x) {
        set_chamber_walls (N, b, x);

        for ( int i = 0; i < boxes.size (); i++ ) {
          set_my_boundary (N, b, x, boxes[i].xy (), boxes[i].z ());
        }

      }

      void gauss_siedel (int N, int b, float * x, float * x0, float a, float c, int * m) {
        auto IX = [=] (int i, int j) { return i + (N + 2)*j; };

        for ( int k = 0; k<20; k++ ) {
          for ( int i = 1; i <= N; i++ ) {
            for ( int j = 1; j <= N; j++ ) {
              if ( !m[IX (i, j)] ) {
                x[IX (i, j)] = (x0[IX (i, j)] + a*(x[IX (i - 1, j)] + x[IX (i + 1, j)] + x[IX (i, j - 1)] + x[IX (i, j + 1)])) / c;
              }
            }
          }
          set_boundary (N, b, x);
        }
      }

      void diffusion (int N, int b, float * x, float * x0, float diff, float dt, int * m) {
        float a = dt * diff * (N * N);
        gauss_siedel (N, b, x, x0, a, 1 + 4 * a, m);
      }

      void advection_step (int N, int b, float * d, float * d0, float * u, float * v, float dt, int * m) {
        auto IX = [=] (int i, int j) { return i + (N + 2)*j; };

        float dt0 = dt*N;
        for ( int i = 1; i <= N; i++ ) {
          for ( int j = 1; j <= N; j++ ) {
            if ( m[IX (i, j)] ) {
              continue;
            }

            // (x, y) is the address to copy from
            float x = i - dt0*u[IX (i, j)], y = j - dt0*v[IX (i, j)];

            // clamp x and y
            if ( x<0.5f ) x = 0.5f; else if ( x>N + 0.5f ) x = N + 0.5f;
            if ( y<0.5f ) y = 0.5f; else if ( y>N + 0.5f ) y = N + 0.5f;

            // s1 and s0 are lerp coordinates [0,1) within the source cell
            int i0 = (int) x, i1 = i0 + 1;
            int j0 = (int) y, j1 = j0 + 1;
            float s1 = x - i0, s0 = 1 - s1;
            float t1 = y - j0, t0 = 1 - t1;

            // sample the source
            d[IX (i, j)] =
              s0*(t0*d0[IX (i0, j0)] + t1*d0[IX (i0, j1)]) +
              s1*(t0*d0[IX (i1, j0)] + t1*d0[IX (i1, j1)])
              ;
          }
        }

        // copy values out to the boundary.
        set_boundary (N, b, d);
      }

      void project (int N, float * u, float * v, float * p, float * div, int * m) {
        auto IX = [=] (int i, int j) { return i + (N + 2)*j; };

        // calculate divergence into div
        // set initial value of p
        for ( int i = 1; i <= N; i++ ) {
          for ( int j = 1; j <= N; j++ ) {
            div[IX (i, j)] = -0.5f*(u[IX (i + 1, j)] - u[IX (i - 1, j)] + v[IX (i, j + 1)] - v[IX (i, j - 1)]) / N;
            p[IX (i, j)] = 0;
          }
        }

        // copy pixels to boundary
        set_boundary (N, 0, div);
        set_boundary (N, 0, p);

        // p += div[x+/-1, y+/-1] * 4;
        gauss_siedel (N, 0, p, div, 1, 4, m);

        // calculate velocity from pressure-like "p"
        for ( int i = 1; i <= N; i++ ) {
          for ( int j = 1; j <= N; j++ ) {
            // u from left and right
            u[IX (i, j)] -= 0.5f*N*(p[IX (i + 1, j)] - p[IX (i - 1, j)]);

            // v from up and down.
            v[IX (i, j)] -= 0.5f*N*(p[IX (i, j + 1)] - p[IX (i, j - 1)]);
          }
        }

        // copy velocity to boundary
        set_boundary (N, 1, u);
        set_boundary (N, 2, v);
      }

      void density_step (int N, float * x, float * x0, float * u, float * v, float diff, float dt, int * m) {
        // apply diffusion to density. If there is no velocity, the value will still spread.
        std::swap (x0, x);
        diffusion (N, 0, x, x0, diff, dt, m);

        // use the velocity field to carry density around.
        std::swap (x0, x);
        advection_step (N, 0, x, x0, u, v, dt, m);
      }

      void velocity_step (int N, float * u, float * v, float * u0, float * v0, float visc, float dt, int * m) {
        // diffuse into neighouring cells
        std::swap (u0, u);
        diffusion (N, 1, u, u0, visc, dt, m);
        std::swap (v0, v);
        diffusion (N, 2, v, v0, visc, dt, m);

        // stabilise the system using poisson
        project (N, u, v, u0, v0, m);
        std::swap (u0, u);
        std::swap (v0, v);

        // use advection to move the velocity itself
        advection_step (N, 1, u, u0, u0, v0, dt, m);
        advection_step (N, 2, v, v0, u0, v0, dt, m);

        // stabilise the system using poisson
        project (N, u, v, u0, v0, m);
      }

      void update (int frame_number) {
        float dt = 1.0f / 30;
        int N = dim.x () - 1;
        assert (density.size () == (N + 2)*(N + 2));
        float *u = vx.data (), *v = vy.data (), *u_prev = prev_vx.data (), *v_prev = prev_vy.data ();
        float *dens = density.data (), *dens_prev = prev_density.data ();
        int* m = mask.data ();
        float visc = 0;
        float diff = 0.0001f;

        //printf("dtot=%f\n", std::accumulate(density.cbegin(), density.cend(), 0.0f));

        // fill the input values
        std::fill (prev_vx.begin (), prev_vx.end (), 0.0f);
        std::fill (prev_vy.begin (), prev_vy.end (), 0.0f);
        std::fill (prev_density.begin (), prev_density.end (), 0.0f);

        // you could use a UI to do this.
        float c = 0;//math::cos(frame_number*0.01f);
        float s = 1;//math::sin (frame_number*0.01f);
        density[fountain_x + (dim.x () + 1) *fountain_y] += 100 * dt;
        u[fountain_x + (dim.x () + 1) *fountain_y] += c * (50 * dt);
        v[fountain_x + (dim.x () + 1) * fountain_y] += s * (50 * dt);

        // step the simulation.
        //get_from_UI( dens_prev, u_prev, v_prev );
        long long t0 = __rdtsc ();
        velocity_step (N, u, v, u_prev, v_prev, visc, dt, m);
        density_step (N, dens, dens_prev, u, v, diff, dt, m);
        long long t1 = __rdtsc ();
        printf ("%lld clocks\n", t1 - t0);

        //printf("dtot=%f\n", std::accumulate(density.cbegin(), density.cend(), 0.0f));

        vertices.resize ((dim.x () + 1)*(dim.y () + 1));
        int stride = (dim.x () + 1);
        size_t d = 0;
        for ( int i = 0; i <= dim.x (); ++i ) {
          for ( int j = 0; j <= dim.y (); ++j ) {
            my_vertex v;
            v.pos = vec3p (i * sx + cx, j * sy + cy, 0);
            float color_value = std::max (0.0f, std::min (density[i + j*stride], 1.0f));
            v.color = vec3p (color_value / 3.0f, color_value, color_value / 3.0f);
            vertices[d++] = v;
          }
        }

        mesh::set_vertices<my_vertex> (vertices);
      }
    };

    ref<mesh_fluid> the_mesh;
  public:
    my_chamber (int argc, char **argv) : app (argc, argv) {}

    void app_init () {
      //init scene
      app_scene = new visual_scene ();
      app_scene->create_default_camera_and_lights ();

      //init fluid_mesh
      material *green = new material (vec4 (1, 0, 0, 1), new param_shader ("shaders/simple_color.vs", "shaders/simple_color.fs"));
      the_mesh = new mesh_fluid (aabb (vec3 (0), vec3 (15)), ivec3 (grid_size, grid_size, 0));
      scene_node *node = new scene_node ();
      app_scene->add_child (node);
      app_scene->add_mesh_instance (new mesh_instance (node, the_mesh, green));

      //init box_mesh
      mat4t sprite_transform;
      sprite_transform.loadIdentity ();
      sprite_transform.translate (vec3 (0, 0, 6));

      image *img = new image ("assets/projects/my_chamber/box.gif");
      material *box_mat = new material (img);
      ref<mesh_sprite> box = new mesh_sprite (vec3 (4, 4, 0), vec2 (7, 7), sprite_transform);
      ref<mesh_sprite> box2 = new mesh_sprite (vec3 (4, 2, 0), vec2 (2, 2), sprite_transform);
      ref<mesh_sprite> box3 = new mesh_sprite (vec3 (4, 4, 0), vec2 (2, 2), sprite_transform);
      //ref<mesh_box> box = new mesh_box(vec3 (2, 2, 2), sprite_transform);
      node = new scene_node ();
      app_scene->add_child (node);
      //app_scene->add_mesh_instance (new mesh_instance (node, box, box_mat));  
      //todo: remove this dumb N argument
      the_mesh->add_box (99, 20, 70, 10);
      the_mesh->add_box (99, 30, 70, 10);
      the_mesh->add_box (99, 40, 70, 10);
      the_mesh->add_box (99, 50, 70, 10);
      the_mesh->add_box (99, 60, 70, 10);
    }

    void draw_world (int x, int y, int w, int h) {
      //todo: check what happens when you add boundaries to the mask
      int vx = 0, vy = 0;
      get_viewport_size (vx, vy);
      app_scene->begin_render (vx, vy, vec4 (0, 0, 0, 1));

      the_mesh->update (get_frame_number ());

      app_scene->update (1.0f / 30);

      app_scene->render ((float) vx / vy);
      
      if(is_key_down(key_up)  ){
        the_mesh->move_box (99, 0, 0, 1);
      }
      if ( is_key_down (key_down) ) {
        the_mesh->move_box (99, 0, 0, -1);
      }
      if ( is_key_down (key_left) ) {
        the_mesh->move_box (99, 0, -1, 0);
      }
      if ( is_key_down (key_right) ) {
        the_mesh->move_box (99, 0, +1, 0);
      }
    }
  };
}
