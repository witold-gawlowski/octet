namespace octet {
  namespace scene {
    class mesh_sprite : public mesh {
      mat4t transform;

      void init (aabb_in size, mat4t_in transform = mat4t ()) {
        this->transform = transform;
        //todo: remove normals
        set_default_attributes (); 
        set_aabb (size);
        update ();
      }

    public:
      RESOURCE_META (mesh_sprite)
        mesh_sprite () {}

      //todo: change vec3 to vec2
      mesh_sprite (vec3_in size, mat4t_in transform = mat4t ()) {
        init (aabb (vec3 (0, 0, 0), size), transform);
      }

      void set_size (vec3_in size, mat4t_in transform = mat4t ()) {
        init (aabb (vec3 (0, 0, 0), size), transform);
        update ();
      }

      /// Generate mesh from parameters.
      virtual void update () {
        aabb aabb_ = get_aabb ();
        mesh::set_shape<math::aabb, mesh::vertex> (aabb_, transform, 1);
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
*     I should add add this later when preparing game physics. 
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

