//! Pile de 8 Mo pour le fil principal sous Windows (1 à 2 Mo par défaut) :
//! l'initialisation de certains pilotes graphiques en consomme beaucoup.

fn main() {
    let target = std::env::var("TARGET").unwrap_or_default();
    if target.contains("windows") {
        if target.ends_with("-gnu") || target.ends_with("-gnullvm") {
            println!("cargo:rustc-link-arg-bins=-Wl,--stack,8388608");
        } else if target.ends_with("-msvc") {
            println!("cargo:rustc-link-arg-bins=/STACK:8388608");
        }
    }
}
