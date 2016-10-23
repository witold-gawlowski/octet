#include <valarray>
#include <memory>

namespace octet {
  class my_chamber : public app {
    class sprite {
      mat4t modelToWorld;
      float halfWidth;
      float halfHeight;
      int texture;
      bool enabled;
    public:
      sprite() {
        texture = 0;
        enabled = true;
      }

      void init(int _texture, float x, float y, float w, float h) {
        modelToWorld.loadIdentity();
        modelToWorld.translate(x, y, 0);
        halfWidth = w * 0.5f;
        halfHeight = h * 0.5f;
        texture = _texture;
        enabled = true;
      }

      void render(texture_shader &shader, mat4t &cameraToWorld) {
        // invisible sprite... used for gameplay.
        if (!texture) return;

        // build a projection matrix: model -> world -> camera -> projection
        // the projection space is the cube -1 <= x/w, y/w, z/w <= 1
        mat4t modelToProjection = mat4t::build_projection_matrix(modelToWorld, cameraToWorld);

        // set up opengl to draw textured triangles using sampler 0 (GL_TEXTURE0)
        glActiveTexture(GL_TEXTURE0);
        //Calling glBindTexture with target set to GL_TEXTURE_2D or GL_TEXTURE_CUBE_MAP and
        //texture set to the name of the new texture binds the texture name to the target of the current active texture unit.
        glBindTexture(GL_TEXTURE_2D, texture);

        // use "old skool" rendering
        //glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST_MIPMAP_NEAREST);
        //glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
        shader.render(modelToProjection, 0);

        // this is an array of the positions of the corners of the sprite in 3D
        // a straight "float" here means this array is being generated here at runtime.
        float vertices[] = {
          -halfWidth, -halfHeight, 0,
          halfWidth, -halfHeight, 0,
          halfWidth,  halfHeight, 0,
          -halfWidth,  halfHeight, 0,
        };

        // attribute_pos (=0) is position of each corner
        // each corner has 3 floats (x, y, z)
        // there is no gap between the 3 floats and hence the stride is 3*sizeof(float)
        glVertexAttribPointer(attribute_pos, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)vertices);
        glEnableVertexAttribArray(attribute_pos);

        // this is an array of the positions of the corners of the texture in 2D
        static const float uvs[] = {
          0,  0,
          1,  0,
          1,  1,
          0,  1,
        };

        // attribute_uv is position in the texture of each corner
        // each corner (vertex) has 2 floats (x, y)
        // there is no gap between the 2 floats and hence the stride is 2*sizeof(float)
        glVertexAttribPointer(attribute_uv, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(float), (void*)uvs);
        glEnableVertexAttribArray(attribute_uv);

        // finally, draw the sprite (4 vertices)
        glDrawArrays(GL_TRIANGLE_FAN, 0, 4);
      }

      // move the object
      void translate(float x, float y) {
        modelToWorld.translate(x, y, 0);
      }

      // position the object relative to another.
      void set_relative(sprite &rhs, float x, float y) {
        modelToWorld = rhs.modelToWorld;
        modelToWorld.translate(x, y, 0);
      }

      // return true if this sprite collides with another.
      // note the "const"s which say we do not modify either sprite
      bool collides_with(const sprite &rhs) const {
        float dx = rhs.modelToWorld[3][0] - modelToWorld[3][0];
        float dy = rhs.modelToWorld[3][1] - modelToWorld[3][1];
        return
          (fabsf(dx) < halfWidth + rhs.halfWidth) &&
          (fabsf(dy) < halfHeight + rhs.halfHeight)
          ;
      }

      bool is_above(const sprite &rhs, float margin) const {
        float dx = rhs.modelToWorld[3][0] - modelToWorld[3][0];

        return
          (fabsf(dx) < halfWidth + margin);
      }

      bool &is_enabled() {
        return enabled;
      }
    };
    class mesh_fluid : public mesh {
      struct my_vertex {
        vec3p pos;
        vec3p color;
      };

      dynarray<my_vertex> vertices;

      std::vector<float> prev_density;
      std::vector<float> prev_vx;
      std::vector<float> prev_vy;

      std::vector<float> density;
      std::vector<float> vx;
      std::vector<float> vy;

      ivec3 dim;
    public:
      mesh_fluid(aabb_in bb, ivec3_in dim) : mesh(), dim(dim) {
        mesh::set_aabb(bb);

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

      /// clamp edges to same (or negative) value of inner pixels
      /// to act as a barrier
      void set_boundary( int N, int b, float * x ) {
        auto IX = [=](int i, int j) { return i +(N+2)*j; };

        //W: i need to keep distance from boxes so that moving them does not consume fluid
	      for ( int i=1 ; i<=N ; i++ ) {
		      x[IX(0  ,i)] = b==1 ? -x[IX(1,i)] : x[IX(1,i)];
		      x[IX(N+1,i)] = b==1 ? -x[IX(N,i)] : x[IX(N,i)];
		      x[IX(i,0  )] = b==2 ? -x[IX(i,1)] : x[IX(i,1)];
		      x[IX(i,N+1)] = b==2 ? -x[IX(i,N)] : x[IX(i,N)];
	      }

	      x[IX(0  ,0  )] = 0.5f*(x[IX(1,0  )]+x[IX(0  ,1)]);
	      x[IX(0  ,N+1)] = 0.5f*(x[IX(1,N+1)]+x[IX(0  ,N)]);
	      x[IX(N+1,0  )] = 0.5f*(x[IX(N,0  )]+x[IX(N+1,1)]);
	      x[IX(N+1,N+1)] = 0.5f*(x[IX(N,N+1)]+x[IX(N+1,N)]);
      }

      /// solve diffusion equations (propagation into adjacent cells)
      /// by repeated aplication of a weighted average.
      /// Use a fixed number of iterations. 
      /// at the end x should not change (you should test this yourself)
      void gauss_siedel( int N, int b, float * x, float * x0, float a, float c ) {
        auto IX = [=](int i, int j) { return i +(N+2)*j; };

	      for ( int k=0 ; k<20 ; k++ ) {
		      for ( int i=1 ; i<=N ; i++ ) {  
            for ( int j=1 ; j<=N ; j++ ) {
			        x[IX(i,j)] = (x0[IX(i,j)] + a*(x[IX(i-1,j)]+x[IX(i+1,j)]+x[IX(i,j-1)]+x[IX(i,j+1)]))/c;
            }
          }
		      set_boundary( N, b, x );
	      }
      }

      /// calculate diffision by approximating with a weighted average
      void diffusion( int N, int b, float * x, float * x0, float diff, float dt ) {
	      float a = dt * diff * (N * N);
	      gauss_siedel( N, b, x, x0, a, 1+4*a );
      }

      /// carry a quantity (velocity or density) from one cell to another.
      void advection_step( int N, int b, float * d, float * d0, float * u, float * v, float dt ) {
        auto IX = [=](int i, int j) { return i +(N+2)*j; };

	      float dt0 = dt*N;
		    for ( int i=1 ; i<=N ; i++ ) {  
          for ( int j=1 ; j<=N ; j++ ) {
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

      // stablisation step. adjust the velocity to prevent increase in energy
      // in the system.
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

      /// Given a velocity field, carry a value around the simulation
      /// and diffuse the value.
      void density_step( int N, float * x, float * x0, float * u, float * v, float diff, float dt ) {
        // apply diffusion to density. If there is no velocity, the value will still spread.
	      std::swap( x0, x );
        diffusion( N, 0, x, x0, diff, dt );

        // use the velocity field to carry density around.
	      std::swap( x0, x );
        advection_step( N, 0, x, x0, u, v, dt );
      }

      /// Compute the new velocity field.
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
        float visc = 0.0f;
        float diff = 0.0f;

        //printf("dtot=%f\n", std::accumulate(density.cbegin(), density.cend(), 0.0f));

        // fill the input values
        std::fill(prev_vx.begin(), prev_vx.end(), 0.0f);
        std::fill(prev_vy.begin(), prev_vy.end(), 0.0f);
        std::fill(prev_density.begin(), prev_density.end(), 0.0f);

        // you could use a UI to do this.
        float c = math::cos(frame_number*0.01f);
        float s = math::sin(frame_number*0.01f);
        density[50 +(dim.x()+1) * 50] += 100 * dt;
        u[50 +(dim.x()+1) * 50] += c * (100 * dt);
        v[50 +(dim.x()+1) * 50] += s * (100 * dt);

        // step the simulation.
	      //get_from_UI( dens_prev, u_prev, v_prev );
        long long t0 = __rdtsc();
	      velocity_step( N, u, v, u_prev, v_prev, visc, dt );
	      density_step( N, dens, dens_prev, u, v, diff, dt );
        long long t1 = __rdtsc();
        printf("%lld clocks\n", t1-t0);

        //printf("dtot=%f\n", std::accumulate(density.cbegin(), density.cend(), 0.0f));

        aabb bb = mesh::get_aabb();
        float sx = bb.get_half_extent().x()*(2.0f/dim.x());
        float sy = bb.get_half_extent().y()*(2.0f/dim.y());
        float cx = bb.get_center().x() - bb.get_half_extent().x();
        float cy = bb.get_center().y() - bb.get_half_extent().y();
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
    enum {
      hero_sprite = 0,
      first_box_sprite,
      num_sprites
    };
    mat4t cameraToWorld;
    texture_shader texture_shader_;
    ref<mesh_fluid> the_mesh;
    ref<visual_scene> app_scene;
    sprite sprites[num_sprites]; 
  public:
    my_chamber(int argc, char **argv) : app(argc, argv) {}
    void app_init() {
      app_scene =  new visual_scene();
      app_scene->create_default_camera_and_lights();

      cameraToWorld.loadIdentity();
      cameraToWorld.translate(vec3(0, 0, 3));
      texture_shader_.init();

      material *red = new material(vec4(1, 0, 0, 1), new param_shader("shaders/simple_color.vs", "shaders/simple_color.fs"));
      the_mesh = new mesh_fluid(aabb(vec3(0), vec3(15)), ivec3(100, 100, 0));
      scene_node *node = new scene_node();
      app_scene->add_child(node);
      app_scene->add_mesh_instance(new mesh_instance(node, the_mesh, red));

      GLuint hero = resource_dict::get_texture_handle(GL_RGBA, "assets/invaderers/ship.gif");
      sprites[hero_sprite].init(hero, 0, 0, 1, 1);
    }

    void draw_world(int x, int y, int w, int h) {
      int vx = 0, vy = 0;
      //W: why does this happen here? Why cant the viewport be size be passed inside begin_render?
      get_viewport_size(vx, vy);
      app_scene->begin_render(vx, vy, vec4(0, 0, 0, 1));
      //W: correct ++i operator in my_bridge

      //W: Strange: the_mesh update actually hides drawn sprites from visibilty(precisely set_vertices function). Why?
      //the_mesh->update(get_frame_number());
      app_scene->update(1.0f/30);
      //app_scene->render((float)vx / vy);

      for (int i = 0; i < num_sprites; ++i) {
        sprites[i].render(texture_shader_, cameraToWorld);
        
      }
    }
  };
}
