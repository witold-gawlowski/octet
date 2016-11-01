#include <memory>
namespace octet {
  class my_chamber : public app {  
    ref<visual_scene> app_scene;
    static const int grid_size = 40;
    static const int fountain_x = 5;
    static const int fountain_y = 5;
    class mesh_fluid : public mesh {
      struct my_vertex {
        vec3p pos;
        vec3p color;
      };

      //boxes can only be square
      dynarray<ref<mesh_sprite>> boxesp;
      dynarray<my_vertex> vertices;
      std::vector<float> prev_density;
      std::vector<float> prev_vx;
      std::vector<float> prev_vy;
      std::vector<float> density;
      std::vector<float> vx;
      std::vector<float> vy;
      ivec3 dim;
      float cx, cy, sx, sy;
      
    public:
      mesh_fluid(aabb_in bb, ivec3_in dim) : mesh(), dim(dim) {
        mesh::set_aabb(bb);
        sx = bb.get_half_extent ().x ()*(2.0f / dim.x ());
        sy = bb.get_half_extent ().y ()*(2.0f / dim.y ());
        cx = bb.get_center ().x () - bb.get_half_extent ().x ();
        cy = bb.get_center ().y () - bb.get_half_extent ().y ();

        density.resize((dim.x()+1)*(dim.y()+1));
        vx.resize((dim.x()+1)*(dim.y()+1));
        vy.resize((dim.x()+1)*(dim.y()+1));

        prev_density.resize((dim.x()+1)*(dim.y()+1));
        prev_vx.resize((dim.x()+1)*(dim.y()+1));
        prev_vy.resize((dim.x()+1)*(dim.y()+1));

        //density[50 +(dim.x()+1) * 50] = 1;
        //prev_vx[50 +(dim.x()+1) * 50] = 1;

        dynarray<uint32_t> indices;
        int stride = dim.x() + 1;
        for (int i = 0; i < dim.x(); ++i) {
          for (int j = 0; j < dim.y(); ++j) {
            indices.push_back((i+1) +(j+0)*stride);
            indices.push_back((i+0) +(j+1)*stride);
            indices.push_back((i+1) +(j+1)*stride);
            indices.push_back((i+1) +(j+0)*stride);
            indices.push_back((i+0) +(j+0)*stride);
            indices.push_back((i+0) +(j+1)*stride);
          }
        }
        set_indices(indices);
        clear_attributes();
        add_attribute(attribute_pos, 3, GL_FLOAT, 0);
        add_attribute(attribute_color, 3, GL_FLOAT, 12);
      }

      void add_box(ref<mesh_sprite> b){
        boxesp.push_back (b);
      }



       bool inside_box(const int& i, const int& j, const vec2& position, const int& size){
        if ( i<position.x () || j<position.y () || i>position.x () + size || j>position.y () + size )
          return false;
        return true;
      }

      bool inside_any_box(const int &i, const int &j){
        for(int k=0; k<boxesp.size(); ++k){
          vec2 pos = boxesp[k]->get_position ().xy();
          if(inside_box(i, j, vec2((pos.x()-cx)/sx, (pos.y()-cy)/sy), boxesp[k]->get_size().x()/sx)){
            return true;
          }
        }
        return false;
      }

      void set_my_boundary(int N, int b, float * x, const vec2& pos, const int& size){
        auto IX = [=] (int i, int j) { return i + (N + 2)*j; };
        for(int i=0; i<size; i++){
          x[IX (pos.x()-1, pos.y()+i)] = b==1 ? -x[IX (pos.x()-2,  pos.y()+i)] : x[IX (pos.x () - 2, pos.y () + i)];
          x[IX (pos.x()+size, pos.y()+i)] = b==1 ? -x[IX (pos.x()+size+1, pos.y()+i)] : x[IX (pos.x () + size + 1, pos.y () + i)];
          x[IX (pos.x () + i, pos.y ()-1)] = b == 2 ? -x[IX (pos.x () + i, pos.y ()-2)] : x[IX (pos.x () + i, pos.y ()-2)];
          x[IX (pos.x () + i, pos.y () + size)] = b == 2 ? -x[IX (pos.x () + i, pos.y () + size +1)] : x[IX (pos.x () + i, pos.y () + size +1)];
        }
        x[IX (pos.x () - 1, pos.y () - 1) ] = 0.5f*(x[IX (pos.x (), pos.y () - 1)] + x[IX (pos.x () - 1, pos.y ())]);
        x[IX (pos.x () + size, pos.y () - 1) ] = 0.5f*(x[IX (pos.x ()+size, pos.y () - 2)] + x[IX (pos.x ()+size + 1, pos.y ()-1)]);
        x[IX (pos.x () +size, pos.y ()+ size) ] = 0.5f*(x[IX (pos.x ()+size+1, pos.y () +size)] + x[IX (pos.x () +size, pos.y ()+size+1)]);
        x[IX (pos.x ()-1, pos.y () +size) ] = 0.5f*(x[IX (pos.x ()-2, pos.y ()+size)] + x[IX (pos.x () - 1, pos.y ()+size+1)]);
      }

      void set_chamber_walls(int N, int b, float * x)
      {
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
      void set_boundary( int N, int b, float * x ) {
        set_chamber_walls (N, b, x);
        
       for ( int i = 0; i < boxesp.size (); ++i ) {
          vec2 pos = boxesp[i]->get_position ().xy ();
          int size = boxesp[i]->get_size ().x ();
          set_my_boundary (N, b, x, vec2 ((pos.x () - cx ) / sx , (pos.y () - cy) / sy), size);
        }
      }

      void gauss_siedel( int N, int b, float * x, float * x0, float a, float c ) {
        auto IX = [=](int i, int j) { return i +(N+2)*j; };

	      for ( int k=0 ; k<20 ; k++ ) {
		      for ( int i=1 ; i<=N ; i++ ) {  
            for ( int j=1 ; j<=N ; j++ ) {
              if (!inside_any_box (i, j)) {
                x[IX (i, j)] = (x0[IX (i, j)] + a*(x[IX (i - 1, j)] + x[IX (i + 1, j)] + x[IX (i, j - 1)] + x[IX (i, j + 1)])) / c;
              }
            }
          }
		      set_boundary( N, b, x );
	      }
      }

      void diffusion( int N, int b, float * x, float * x0, float diff, float dt ) {
	      float a = dt * diff * (N * N);
	      gauss_siedel( N, b, x, x0, a, 1+4*a );
      }

      void advection_step( int N, int b, float * d, float * d0, float * u, float * v, float dt ) {
        auto IX = [=](int i, int j) { return i +(N+2)*j; };

	      float dt0 = dt*N;
		    for ( int i=1 ; i<=N ; i++ ) {  
          for ( int j=1 ; j<=N ; j++ ) {
            if (inside_any_box (i, j) ) {
              continue;
            }

            // (x, y) is the address to copy from
		        float x = i-dt0*u[IX(i,j)], y = j-dt0*v[IX(i,j)];

            // clamp x and y
		        if (x<0.5f) x=0.5f; else if (x>N+0.5f) x=N+0.5f;
		        if (y<0.5f) y=0.5f; else if (y>N+0.5f) y=N+0.5f;

            // s1 and s0 are lerp coordinates [0,1) within the source cell
            int i0=(int)x, i1=i0+1;
            int j0=(int)y, j1=j0+1;
		        float s1 = x - i0, s0 = 1 - s1;
            float t1 = y - j0, t0 = 1 - t1;

            // sample the source
		        d[IX(i,j)] =
              s0*(t0*d0[IX(i0,j0)]+t1*d0[IX(i0,j1)]) +
              s1*(t0*d0[IX(i1,j0)]+t1*d0[IX(i1,j1)])
            ;
          }
        }

        // copy values out to the boundary.
	      set_boundary( N, b, d );
      }

      void project( int N, float * u, float * v, float * p, float * div ) {
        auto IX = [=](int i, int j) { return i +(N+2)*j; };

        // calculate divergence into div
        // set initial value of p
		    for ( int i=1 ; i<=N ; i++ ) {  
          for ( int j=1 ; j<=N ; j++ ) {
		        div[IX(i,j)] = -0.5f*(u[IX(i+1,j)]-u[IX(i-1,j)]+v[IX(i,j+1)]-v[IX(i,j-1)])/N;
		        p[IX(i,j)] = 0;
          }
        }

        // copy pixels to boundary
	      set_boundary( N, 0, div );
        set_boundary( N, 0, p );

        // p += div[x+/-1, y+/-1] * 4;
	      gauss_siedel( N, 0, p, div, 1, 4 );

        // calculate velocity from pressure-like "p"
		    for ( int i=1 ; i<=N ; i++ ) {  
          for ( int j=1 ; j<=N ; j++ ) {
            // u from left and right
		        u[IX(i,j)] -= 0.5f*N*(p[IX(i+1,j)]-p[IX(i-1,j)]);

            // v from up and down.
		        v[IX(i,j)] -= 0.5f*N*(p[IX(i,j+1)]-p[IX(i,j-1)]);
          }
        }

        // copy velocity to boundary
	      set_boundary( N, 1, u );
        set_boundary( N, 2, v );
      }

      void density_step( int N, float * x, float * x0, float * u, float * v, float diff, float dt ) {
        // apply diffusion to density. If there is no velocity, the value will still spread.
	      std::swap( x0, x );
        diffusion( N, 0, x, x0, diff, dt );

        // use the velocity field to carry density around.
	      std::swap( x0, x );
        advection_step( N, 0, x, x0, u, v, dt );
      }

      void velocity_step( int N, float * u, float * v, float * u0, float * v0, float visc, float dt ) {
        // diffuse into neighouring cells
	      std::swap( u0, u );
        diffusion( N, 1, u, u0, visc, dt );
	      std::swap( v0, v );
        diffusion( N, 2, v, v0, visc, dt );

        // stabilise the system using poisson
	      project( N, u, v, u0, v0 );
	      std::swap( u0, u );
        std::swap( v0, v );

        // use advection to move the velocity itself
	      advection_step( N, 1, u, u0, u0, v0, dt );
        advection_step( N, 2, v, v0, u0, v0, dt );

        // stabilise the system using poisson
	      project( N, u, v, u0, v0 );
      }

      void update(int frame_number) {
        float dt = 1.0f / 30;
        int N = dim.x()-1;
        assert(density.size() == (N+2)*(N+2));
        float *u = vx.data(), *v = vy.data(), *u_prev = prev_vx.data(), *v_prev = prev_vy.data();
        float *dens = density.data(), *dens_prev = prev_density.data();
        float visc = 0;
        float diff = 0.0000f;

        //printf("dtot=%f\n", std::accumulate(density.cbegin(), density.cend(), 0.0f));

        // fill the input values
        std::fill(prev_vx.begin(), prev_vx.end(), 0.0f);
        std::fill(prev_vy.begin(), prev_vy.end(), 0.0f);
        std::fill(prev_density.begin(), prev_density.end(), 0.0f);

        // you could use a UI to do this.
        float c = 1;//math::cos(frame_number*0.01f);
        float s = 1;// math::sin (frame_number*0.01f);
        density[fountain_x +(dim.x()+1) *fountain_y] += 100 * dt;
        u[fountain_x +(dim.x()+1) *fountain_y] += c * (50 * dt);
        v[fountain_x + (dim.x()+1) * fountain_y] += s * (50 * dt);

        // step the simulation.
	      //get_from_UI( dens_prev, u_prev, v_prev );
        long long t0 = __rdtsc();
	      velocity_step( N, u, v, u_prev, v_prev, visc, dt );
	      density_step( N, dens, dens_prev, u, v, diff, dt );
        long long t1 = __rdtsc();
        printf("%lld clocks\n", t1-t0);

        //printf("dtot=%f\n", std::accumulate(density.cbegin(), density.cend(), 0.0f));

        vertices.resize((dim.x()+1)*(dim.y()+1));
        int stride =(dim.x()+1);
        size_t d = 0;
        for (int i = 0; i <= dim.x(); ++i) {
          for (int j = 0; j <= dim.y(); ++j) {
            my_vertex v;
            v.pos = vec3p(i * sx + cx, j * sy + cy, 0);
            float color_value = std::max(0.0f, std::min(density[i + j*stride], 1.0f));
            v.color = vec3p(color_value/3.0f, color_value, color_value/3.0f);
            vertices[d++] = v;
          }
        }

        mesh::set_vertices<my_vertex>(vertices);
      }
    };

    ref<mesh_fluid> the_mesh;
  public:
    my_chamber(int argc, char **argv) : app(argc, argv) {
    }

    void app_init() {
      //init scene
      app_scene =  new visual_scene();
      app_scene->create_default_camera_and_lights();
      
      //init fluid_mesh
      material *green = new material(vec4(1, 0, 0, 1), new param_shader("shaders/simple_color.vs", "shaders/simple_color.fs"));
      the_mesh = new mesh_fluid(aabb(vec3(0), vec3(15)), ivec3(grid_size, grid_size, 0));
      scene_node *node = new scene_node();
      app_scene->add_child(node);
      app_scene->add_mesh_instance(new mesh_instance(node, the_mesh, green));

      //init box_mesh
      mat4t sprite_transform;
      sprite_transform.loadIdentity ();
      sprite_transform.translate (vec3 (0, 0, 6));

      image *img = new image ("assets/projects/my_chamber/box.gif");
      material *box_mat= new material (img);
      ref<mesh_sprite> box = new mesh_sprite(vec3 (4, 4, 0), vec2(7, 7), sprite_transform);
      ref<mesh_sprite> box2 = new mesh_sprite (vec3 (4, 2, 0), vec2 (2, 2), sprite_transform);
      ref<mesh_sprite> box3 = new mesh_sprite (vec3 (4, 4, 0), vec2 (2, 2), sprite_transform);
      //ref<mesh_box> box = new mesh_box(vec3 (2, 2, 2), sprite_transform);
      node = new scene_node ();
      app_scene->add_child (node);
      app_scene->add_mesh_instance (new mesh_instance (node, box, box_mat));  
      the_mesh->add_box (box);
      //the_mesh->add_box (box2);
      //the_mesh->add_box (box3);
    }

    void draw_world(int x, int y, int w, int h) {
      int vx = 0, vy = 0;
      get_viewport_size(vx, vy);
      app_scene->begin_render(vx, vy, vec4(0, 0, 0, 1));

      the_mesh->update(get_frame_number());

      app_scene->update(1.0f/30);

      app_scene->render((float)vx / vy);
    }
  };
}
