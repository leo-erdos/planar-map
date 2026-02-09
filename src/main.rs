use planar_map::{PlanarMap, labels, tree_to_quadrangulation};
use rand::{SeedableRng, rngs::StdRng};

fn main() {
    // 1, 3: Doesn't like rev
    // 0, 2: Likes rev
    let mut rng = StdRng::from_seed([6u8; 32]);
    let x = PlanarMap::random_tree(10000000, &mut rng);
    let labels = labels(&x, &mut rng);
    let quad = tree_to_quadrangulation(&x, &labels);
    dbg!(quad.faces().first());
}
