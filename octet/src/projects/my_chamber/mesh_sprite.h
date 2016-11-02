namespace octet {
  namespace scene {
    class mesh_sprite : public mesh {
      mat4t transform;
      quad _quad;

    public:
      RESOURCE_META (mesh_sprite)
        mesh_sprite () {}

      /*  todo: override: 
      *   virtual btCollisionShape *get_bullet_shape() - mesh member   
      */ 
      
      vec3 get_position(){ return _quad.get_position () * transform; }

      vec2 get_size () { return _quad.get_half_extent () * 2; }

      vec3 set_position (vec3 p) { _quad.set_position (p); }

      vec2 set_size (vec2 s) { _quad.set_half_extent (s); }

      mesh_sprite (vec3 pos, vec2 size, mat4t transform) {
        set_default_attributes ();
        this->transform = transform;
        _quad = quad (pos, size / 2);
        vec3 corners[4] = {
          vec3(pos + size / 2),
          vec3(pos - size / 2),
          vec3(pos + size.rot90() / 2),
          vec3(pos - size.rot90() / 2)
        };
        set_aabb (aabb (corners, corners+3));
        aabb aabb_ = get_aabb ();
        //printf ("sprite aabb: %f, %f\n", aabb_.get_center(), aabb_.get_half_extent());
        //std::cout << "one of them: " << vec3(size) << std::endl;
        //std::cout << "corners" << std::endl << corners[0] << std::endl << corners[1] << std::endl << corners[2] <<
        //std::endl << corners[3] << std::endl;
        //set_aabb (aabb(pos, vec3(size.x())));
        update ();
      }

      /// Generate mesh from parameters.
      virtual void update () {
        aabb aabb_ = get_aabb ();
        mesh::set_shape<math::quad, mesh::vertex> (_quad, transform, 1);
        //dump(log("zzz\n"));
      }

      /// Serialise the box
      void visit (visitor &v) {
        mesh::visit (v);
        if ( v.is_reader () ) {
          update ();
        }
      }

/*
*     todo: add add this later when preparing game physics. 
*
#ifdef OCTET_BULLET
      /// Get a bullet shape object for this mesh
      btCollisionShape *get_bullet_shape () {
        return new btBoxShape (get_btVector3 (get_aabb ().get_half_extent ()));
      }

      /// Get a static bullet shape object for this mesh
      btCollisionShape *get_static_bullet_shape () {
        return new btBoxShape (get_btVector3 (get_aabb ().get_half_extent ()));
      }
#endif
*/
    };
  }
}

