use bytemuck::cast_slice;
use std::{collections::HashSet, fs::File, hint::black_box, io::Write, time::Instant};

use planar_map::{
    HalfEdgeColor, PlanarMap, augmented_distance_matrix, distance_matrix, initial_3d_embedding,
    labels, quadrangulation_to_triangulation, save_obj, save_quad_obj, schnyder_embedding,
    schnyder_woods, simplify_triangulation, tree_to_quadrangulation,
};
use rand::{SeedableRng, rngs::StdRng};

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let faces: usize = args[1].parse().unwrap();
    let mut rng = StdRng::from_seed([1u8; 32]);
    let mut map = PlanarMap::random_tree(faces, &mut rng);
    let labels = labels(&map, &mut rng);
    tree_to_quadrangulation(&mut map, &labels);
    let face_walk = map.rotation_system().compose(map.edge_bindings());
    for face in face_walk.cycles() {
        assert_eq!(face.count(), 4);
    }
    let quad_dst = distance_matrix(&map);
    let n_vertices_no_augment: u64 = map.rotation_system().cycles().count() as u64;
    let diag_start = map.num_half_edges();
    quadrangulation_to_triangulation(&mut map);
    let n_half_edges_no_augment = map.num_half_edges();
    let true_edges = simplify_triangulation(&mut map);
    let aug_dst = augmented_distance_matrix(
        &map,
        &quad_dst,
        n_vertices_no_augment,
        n_half_edges_no_augment,
    );
    let n_vertices = map.rotation_system().cycles().count() as u64;
    let colors = schnyder_woods(&map);
    let embedding = schnyder_embedding(&map, &colors);
    let vertex_map = map.rotation_system().cycle_map();
    let initial_embedding = initial_3d_embedding(
        &embedding,
        embedding[&0].1 as u64,
        &mut rng,
        &HashSet::from([
            0,
            vertex_map[map.rotation_system().at(map.edge_bindings().at(0)) as usize],
            vertex_map[map.rotation_system().at(map
                .edge_bindings()
                .at(map.rotation_system().at(map.edge_bindings().at(0))))
                as usize],
        ]),
    );
    println!(
        "Faces {}",
        map.rotation_system()
            .compose(map.edge_bindings())
            .cycles()
            .count()
    );
    println!("Vertices {}", map.rotation_system().cycles().count());
    println!("Edges {}", map.edges().count());
    save_obj(&map, &initial_embedding, "small.obj");
    let mut dist_file = File::create("distances").unwrap();
    for d in aug_dst {
        write!(dist_file, "{d}\n").unwrap();
    }
    save_quad_obj(&map, diag_start, &true_edges, "quad_edges.obj");
}
