use crate::Polygon;

pub(crate) mod transform{
    use crate::{Polygon, PolygonOwned};

    pub(crate) fn rotate(polygon: Polygon, angle: f64) -> PolygonOwned {
        polygon
            .iter()
            .map(|&(x, y)| {
                (
                    x * angle.cos() - y * angle.sin(),
                    x * angle.sin() + y * angle.cos(),
                )
            })
            .collect::<Vec<_>>()
    }
    
    pub(crate) fn translate(polygon: Polygon, x: f64, y: f64) -> PolygonOwned {
        polygon
            .iter()
            .map(|&(x0, y0)| (x0 + x, y0 + y))
            .collect::<Vec<_>>()
    }
    
    pub(crate) fn scale(polygon: Polygon, x: f64, y: f64) -> PolygonOwned {
        polygon
            .iter()
            .map(|&(x0, y0)| (x0 * x, y0 * y))
            .collect::<Vec<_>>()
    }
}


pub(crate) fn find_horizontal_intersect(polygon: Polygon, y: f64) -> Vec<f64> {
    let mut intersections = Vec::new();
    for i in 0..polygon.len() {
        let j = (i + 1) % polygon.len();
        let (x1, y1) = polygon[i];
        let (x2, y2) = polygon[j];

        if (y1 - y) * (y2 - y) < 0.0 {
            let x = x1 + (x2 - x1) * (y - y1) / (y2 - y1);
            intersections.push(x);
        }
    }
    intersections
}

pub(crate) fn find_bottom(polygon: Polygon) -> usize {
    let mut min_x = polygon[0].0;
    let mut min_y = polygon[0].1;
    let mut min = 0;
    for i in 1..polygon.len() {
        let (x, y) = polygon[i];
        if y < min_y || (y == min_y && x < min_x) {
            min_x = x;
            min_y = y;
            min = i;
        }
    }

    min
}