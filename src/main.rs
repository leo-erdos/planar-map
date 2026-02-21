use planar_map::{PlanarMap, labels, quadrangulation_to_triangulation, tree_to_quadrangulation};
use rand::{SeedableRng, rngs::StdRng};

fn main() {
    let mut rng = StdRng::from_seed([6u8; 32]);
    let mut map = PlanarMap::random_tree(100, &mut rng);
    let labels = labels(&map, &mut rng);
    tree_to_quadrangulation(&mut map, &labels);
    dbg!(map.faces());
}
