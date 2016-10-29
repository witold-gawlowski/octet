////////////////////////////////////////////////////////////////////////////////
//
// (C) Andy Thomason 2012-2014
//
// Modular Framework for OpenGLES2 rendering on multiple platforms.
//
// Axis aligned bounding box
//

namespace octet {
  namespace math {
    class quad {
      vec3 position;
      vec2 half_extent;
      vec3 axes[3];
    public:
      quad () {
        position = vec3 (0, 0, 0);
        half_extent = vec2 (0.5f, 0.5f);
        axes[0] = vec3 (1, 0, 0);
        axes[1] = vec3 (0, 1, 0);
        axes[2] = vec3 (0, 0, 1);
      }

      quad (const vec3 &center_, const vec2 &half_extent_, const mat4t &transform = mat4t()) {
        position = center_* transform;
        half_extent = half_extent_;
        axes[0] = transform[0].xyz ();
        axes[1] = transform[1].xyz ();
        axes[2] = transform[2].xyz ();
      }

      // Get the center of the bounding box
      const vec3 get_position () const {
        return position;
      }

      const vec2 get_half_extent () const {
        return half_extent;
      }

      //todo: implement later
      
       /*const char *toString (char *dest, size_t len) const {
        char tmp[2][64];
        snprintf (dest, len, "[%s, %s]", position.toString (tmp[0], sizeof (tmp[0])), half_extent.toString (tmp[1], sizeof (tmp[1])));
        return dest;
      }*/




      template <class sink_t> void get_geometry (sink_t &sink, int) {
        static const float vertices[4 * 6 * 8] = {
          // front
#undef OCTET_FACE
#define OCTET_FACE(X, Y, U, V) X, Y, 1,  0, 0, 1, U, V,
          OCTET_FACE (1,  1,  1, 1) OCTET_FACE (1, -1,  1, 0) OCTET_FACE (-1, -1,  0, 0) OCTET_FACE (-1,  1,  0, 1)

          // back
#undef OCTET_FACE
#define OCTET_FACE(X, Y, U, V) -(X), (Y), -1,  0, 0, -1, U, V,
          OCTET_FACE (1,  1,  1, 1) OCTET_FACE (1, -1,  1, 0) OCTET_FACE (-1, -1,  0, 0) OCTET_FACE (-1,  1,  0, 1)

          // top
#undef OCTET_FACE
#define OCTET_FACE(X, Y, U, V) Y, 1, X,  0, 1, 0, U, V,
          OCTET_FACE (1,  1,  1, 1) OCTET_FACE (1, -1,  1, 0) OCTET_FACE (-1, -1,  0, 0) OCTET_FACE (-1,  1,  0, 1)

          // bottom
#undef OCTET_FACE
#define OCTET_FACE(X, Y, U, V) (Y), -1, -(X),  0, 0, -1, U, V,
          OCTET_FACE (1,  1,  1, 1) OCTET_FACE (1, -1,  1, 0) OCTET_FACE (-1, -1,  0, 0) OCTET_FACE (-1,  1,  0, 1)

          // right
#undef OCTET_FACE
#define OCTET_FACE(X, Y, U, V) 1, X, Y,  1, 0, 0, U, V,
          OCTET_FACE (1,  1,  1, 1) OCTET_FACE (1, -1,  1, 0) OCTET_FACE (-1, -1,  0, 0) OCTET_FACE (-1,  1,  0, 1)

          // left
#undef OCTET_FACE
#define OCTET_FACE(X, Y, U, V) -1, -(X), (Y),  -1, 0, 0, U, V,
          OCTET_FACE (1,  1,  1, 1) OCTET_FACE (1, -1,  1, 0) OCTET_FACE (-1, -1,  0, 0) OCTET_FACE (-1,  1,  0, 1)

#undef OCTET_FACE
        };

        // 3 0
        // 2 1
        static const uint32_t indices[6 * 6] = {
          0 + 4 * 0, 1 + 4 * 0, 3 + 4 * 0,  1 + 4 * 0, 2 + 4 * 0, 3 + 4 * 0,
          0 + 4 * 1, 1 + 4 * 1, 3 + 4 * 1,  1 + 4 * 1, 2 + 4 * 1, 3 + 4 * 1,
          0 + 4 * 2, 1 + 4 * 2, 3 + 4 * 2,  1 + 4 * 2, 2 + 4 * 2, 3 + 4 * 2,
          0 + 4 * 3, 1 + 4 * 3, 3 + 4 * 3,  1 + 4 * 3, 2 + 4 * 3, 3 + 4 * 3,
          0 + 4 * 4, 1 + 4 * 4, 3 + 4 * 4,  1 + 4 * 4, 2 + 4 * 4, 3 + 4 * 4,
          0 + 4 * 5, 1 + 4 * 5, 3 + 4 * 5,  1 + 4 * 5, 2 + 4 * 5, 3 + 4 * 5,
        };

        sink.reserve (6 * 4, 6 * 6);

        for ( unsigned i = 0; i != 6 * 4; ++i ) {
          const float *fs = vertices + i * 8;
          sink.add_vertex (position + vec3 (fs[0], fs[1], fs[2]) * half_extent, vec3 (fs[3], fs[4], fs[5]), vec3 (fs[6], fs[7], 0));
        }

        for ( unsigned i = 0; i != 6 * 6; i += 3 ) {
          sink.add_triangle (indices[i + 0], indices[i + 1], indices[i + 2]);
        }
      }
    };

   /* std::ostream &operator <<(std::ostream &os, quad rhs) {
      char tmp[256];
      os << rhs.toString (tmp, sizeof (tmp));
      return os;
    }*/
  }
}

