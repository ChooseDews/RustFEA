use std::borrow::Borrow;

use nalgebra::{self as na, Vector6};
use na::{DVector};
use serde::{Serialize, Deserialize};
use log::{debug, info, warn};
use toml::Value;

pub fn print_max_displacement(displacement: &DVector<f64>) {
    let max_disp = displacement.iter().max_by(|a, b| a.partial_cmp(b).unwrap()).unwrap();
    debug!("Max displacement: {}", max_disp);
}

pub fn check_for_nans<I>(iter: I) -> bool
where
    I: IntoIterator,
    I::Item: std::borrow::Borrow<f64>,
{
    for (i, x) in iter.into_iter().enumerate() {
        if x.borrow().is_nan() {
            debug!("NaN at index: {}", i);
            return true;
        }
    }
    false
}


pub fn safe_component_div(a: &DVector<f64>, b: &DVector<f64>) -> DVector<f64> {
    let mut result = a.clone();
    for (i, &x) in a.iter().enumerate() {
        if b[i] != 0.0 {
            result[i] = x / b[i];
        }else{
            if x > 0.0 {
                warn!("Zero Mass with an applied force at node: {}. What does it mean? zero hopefully.", i);
            }
            result[i] = 0.0;
        }
    }
    result
}


pub fn compute_von_mises(stress_vector: Vector6<f64>) -> f64 {
    //compute and return the von mises stress here
    let s11 = stress_vector[0];
    let s22 = stress_vector[1];
    let s33 = stress_vector[2];
    let s12 = stress_vector[3];
    let s23 = stress_vector[4];
    let s13 = stress_vector[5];
    let s_vm = ((s11 - s22).powi(2) + (s22 - s33).powi(2) + (s33 - s11).powi(2) + 6.0 * (s12.powi(2) + s23.powi(2) + s13.powi(2))) / 2.0;
    s_vm.sqrt()
}

//compute the y displacement of a cantilever beam with a point load at the end
pub fn eular_beam_displacement(x: f64, length: f64, load: f64, modulus: f64, moment_of_inertia: f64) -> f64 {
    let y = load * x.powi(2) * (3.0 * length - x) / (6.0 * modulus * moment_of_inertia);
    y
}

#[derive(Serialize, Deserialize, Clone)]
pub struct Keyword {
    keyword: String,
    pub value: Value
}

impl Keyword {
    pub fn print(&self) {
        debug!("{:}={:?}", self.keyword, self.value);
    }
}

#[derive(Serialize, Deserialize, Clone)]
pub struct Keywords {
    values: Vec<Keyword>,
}

impl Keywords {
    pub fn new() -> Self {
        Keywords {
            values: Vec::new()
        }
    }
    pub fn add_keyword(&mut self, keyword: String, value: Value) {
        self.values.push(Keyword { keyword: keyword.to_uppercase(), value });
    }
    pub fn get_keyword(&self, keyword: &str) -> Option<&Keyword> {
        self.values.iter().find(|k| k.keyword == keyword.to_uppercase())
    }
    pub fn keyword_exists(&self, keyword: &str) -> bool {
        self.values.iter().any(|k| k.keyword == keyword.to_uppercase())
    }   
    pub fn get_keywords(&self, keyword: &str) -> Vec<&Keyword> { //some keywords may have multiple values or instances
        self.values.iter().filter(|k| k.keyword == keyword.to_uppercase()).collect()
    }
    pub fn get_float(&self, keyword: &str) -> Option<f64> {
        if let Some(keyword) = self.values.iter().find(|k| k.keyword == keyword.to_uppercase()) {
            //if int convert to float
            if keyword.value.is_integer() {
                return Some(keyword.value.as_integer().unwrap() as f64);
            }else if keyword.value.is_float() {
                return Some(keyword.value.as_float().unwrap());
            }
        }
        None
    }
    pub fn get_int(&self, keyword: &str) -> Option<i64> {
        if let Some(keyword) = self.values.iter().find(|k| k.keyword == keyword.to_uppercase()) {
            return Some(keyword.value.as_integer().unwrap());
        }
        None
    }
    pub fn get_string(&self, keyword: &str) -> Option<String> {
        if let Some(keyword) = self.values.iter().find(|k| k.keyword == keyword.to_uppercase()) {
            return Some(keyword.value.as_str().unwrap().to_string());
        }
        None
    }
    pub fn get_bool(&self, keyword: &str) -> Option<bool> {
        if let Some(keyword) = self.values.iter().find(|k| k.keyword == keyword.to_uppercase()) {
            return Some(keyword.value.as_bool().unwrap());
        }
        None
    }
    pub fn add(&mut self, keyword: &str, value: Value) {
        let keyword = Keyword { 
            keyword: keyword.to_uppercase(), 
            value: value
        };
        self.values.push(keyword);
    }

    pub fn print(&self) {
        debug!("Keywords ({}):", self.values.len());
        for keyword in &self.values {
            keyword.print();
        }
    }
}