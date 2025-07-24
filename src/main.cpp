#include "AtomicOrbital.hpp"
#include "PrimitiveGaussian.hpp"
#include "Utils.hpp"

#include <iostream>

int main()
{
    // Water molecule (H2O) example with STO-3G basis set
    const PrimitiveGaussian H1_1S_pg1(3.425250914, 0.1543289673, Vec3(0.0, 1.43233673, -0.96104039), 0, 0, 0);
    const PrimitiveGaussian H1_1S_pg2(0.6239137298, 0.5353281423, Vec3(0.0, 1.43233673, -0.96104039), 0, 0, 0);
    const PrimitiveGaussian H1_1S_pg3(0.1688554040, 0.4446345422, Vec3(0.0, 1.43233673, -0.96104039), 0, 0, 0);
    const AtomicOrbital H1_1S(Vec3(0.0, 1.43233673, -0.96104039), {H1_1S_pg1, H1_1S_pg2, H1_1S_pg3});

    const PrimitiveGaussian H2_1S_pg1(3.425250914, 0.1543289673, Vec3(0.0, -1.43233673, -0.96104039), 0, 0, 0);
    const PrimitiveGaussian H2_1S_pg2(0.6239137298, 0.5353281423, Vec3(0.0, -1.43233673, -0.96104039), 0, 0, 0);
    const PrimitiveGaussian H2_1S_pg3(0.1688554040, 0.4446345422, Vec3(0.0, -1.43233673, -0.96104039), 0, 0, 0);
    const AtomicOrbital H2_1S(Vec3(0.0, -1.43233673, -0.96104039), {H2_1S_pg1, H2_1S_pg2, H2_1S_pg3});

    const PrimitiveGaussian O_1S_pg1(130.7093214, 0.1543289673, Vec3(0.0, 0.0, 0.24026010), 0, 0, 0);
    const PrimitiveGaussian O_1S_pg2(23.80886605, 0.5353281423, Vec3(0.0, 0.0, 0.24026010), 0, 0, 0);
    const PrimitiveGaussian O_1S_pg3(6.443608313, 0.4446345422, Vec3(0.0, 0.0, 0.24026010), 0, 0, 0);
    const AtomicOrbital O_1S(Vec3(0.0, 0.0, 0.24026010), {O_1S_pg1, O_1S_pg2, O_1S_pg3});

    const PrimitiveGaussian O_2S_pg1(5.033151319, -0.09996722919, Vec3(0.0, 0.0, 0.24026010), 0, 0, 0);
    const PrimitiveGaussian O_2S_pg2(1.169596125, 0.3995128261, Vec3(0.0, 0.0, 0.24026010), 0, 0, 0);
    const PrimitiveGaussian O_2S_pg3(0.3803889600, 0.7001154689, Vec3(0.0, 0.0, 0.24026010), 0, 0, 0);
    const AtomicOrbital O_2S(Vec3(0.0, 0.0, 0.24026010), {O_2S_pg1, O_2S_pg2, O_2S_pg3});

    const PrimitiveGaussian O_2Px_pg1(5.033151319, 0.1559162750, Vec3(0.0, 0.0, 0.24026010), 1, 0, 0);
    const PrimitiveGaussian O_2Px_pg2(1.169596125, 0.6076837186, Vec3(0.0, 0.0, 0.24026010), 1, 0, 0);
    const PrimitiveGaussian O_2Px_pg3(0.3803889600, 0.3919573931, Vec3(0.0, 0.0, 0.24026010), 1, 0, 0);
    const AtomicOrbital O_2Px(Vec3(0.0, 0.0, 0.24026010), {O_2Px_pg1, O_2Px_pg2, O_2Px_pg3});

    const PrimitiveGaussian O_2Py_pg1(5.033151319, 0.1559162750, Vec3(0.0, 0.0, 0.24026010), 0, 1, 0);
    const PrimitiveGaussian O_2Py_pg2(1.169596125, 0.6076837186, Vec3(0.0, 0.0, 0.24026010), 0, 1, 0);
    const PrimitiveGaussian O_2Py_pg3(0.3803889600, 0.3919573931, Vec3(0.0, 0.0, 0.24026010), 0, 1, 0);
    const AtomicOrbital O_2Py(Vec3(0.0, 0.0, 0.24026010), {O_2Py_pg1, O_2Py_pg2, O_2Py_pg3});

    const PrimitiveGaussian O_2Pz_pg1(5.033151319, 0.1559162750, Vec3(0.0, 0.0, 0.24026010), 0, 0, 1);
    const PrimitiveGaussian O_2Pz_pg2(1.169596125, 0.6076837186, Vec3(0.0, 0.0, 0.24026010), 0, 0, 1);
    const PrimitiveGaussian O_2Pz_pg3(0.3803889600, 0.3919573931, Vec3(0.0, 0.0, 0.24026010), 0, 0, 1);
    const AtomicOrbital O_2Pz(Vec3(0.0, 0.0, 0.24026010), {O_2Pz_pg1, O_2Pz_pg2, O_2Pz_pg3});

    // H2 molecule example with STO-3G basis set
    // const PrimitiveGaussian H1_1S_pg1(3.425250914, 0.1543289673, Vec3(0.0, 0.0, 0.0), 0, 0, 0);
    // const PrimitiveGaussian H1_1S_pg2(0.6239137298, 0.5353281423, Vec3(0.0, 0.0, 0.0), 0, 0, 0);
    // const PrimitiveGaussian H1_1S_pg3(0.1688554040, 0.4446345422, Vec3(0.0, 0.0, 0.0), 0, 0, 0);
    // const AtomicOrbital H1_1S(Vec3(0.0, 0.0, 0.0), {H1_1S_pg1, H1_1S_pg2, H1_1S_pg3});
    //
    // const PrimitiveGaussian H2_1S_pg1(3.425250914, 0.1543289673, Vec3(0, 0, 1.4), 0, 0, 0);
    // const PrimitiveGaussian H2_1S_pg2(0.6239137298, 0.5353281423, Vec3(0, 0, 1.4), 0, 0, 0);
    // const PrimitiveGaussian H2_1S_pg3(0.1688554040, 0.4446345422, Vec3(0, 0, 1.4), 0, 0, 0);
    // const AtomicOrbital H2_1S(Vec3(0.0, 0.0, 1.4), {H2_1S_pg1, H2_1S_pg2, H2_1S_pg3});

    // H2 molecule example with STO-3G basis set
    // const PrimitiveGaussian H_1S_pg1(3.425250914, 0.1543289673, Vec3(0.0, 0.0, -0.1527810236), 0, 0, 0);
    // const PrimitiveGaussian H_1S_pg2(0.6239137298, 0.5353281423, Vec3(0.0, 0.0, -0.1527810236), 0, 0, 0);
    // const PrimitiveGaussian H_1S_pg3(0.1688554040, 0.4446345422, Vec3(0.0, 0.0, -0.1527810236), 0, 0, 0);
    // const AtomicOrbital H_1S(Vec3(0.0, 0.0, -0.1527810236), {H_1S_pg1, H_1S_pg2, H_1S_pg3});
    //
    // const PrimitiveGaussian F_1S_pg1(166.6791340, 0.1543289673, Vec3(0.0, 0.0, 1.6527810236), 0, 0, 0);
    // const PrimitiveGaussian F_1S_pg2(30.36081233, 0.5353281423, Vec3(0.0, 0.0, 1.6527810236), 0, 0, 0);
    // const PrimitiveGaussian F_1S_pg3(8.216820672, 0.4446345422, Vec3(0.0, 0.0, 1.6527810236), 0, 0, 0);
    // const AtomicOrbital F_1S(Vec3(0.0, 0.0, 1.6527810236), {F_1S_pg1, F_1S_pg2, F_1S_pg3});
    //
    // const PrimitiveGaussian F_2S_pg1(6.464803249, -0.09996722919, Vec3(0.0, 0.0, 1.6527810236), 0, 0, 0);
    // const PrimitiveGaussian F_2S_pg2(1.502281245, 0.3995128261, Vec3(0.0, 0.0, 1.6527810236), 0, 0, 0);
    // const PrimitiveGaussian F_2S_pg3(0.4885884864, 0.7001154689, Vec3(0.0, 0.0, 1.6527810236), 0, 0, 0);
    // const AtomicOrbital F_2S(Vec3(0.0, 0.0, 1.6527810236), {F_2S_pg1, F_2S_pg2, F_2S_pg3});
    //
    // const PrimitiveGaussian F_2Px_pg1(6.464803249, 0.1559162750, Vec3(0.0, 0.0, 1.6527810236), 1, 0, 0);
    // const PrimitiveGaussian F_2Px_pg2(1.502281245, 0.6076837186, Vec3(0.0, 0.0, 1.6527810236), 1, 0, 0);
    // const PrimitiveGaussian F_2Px_pg3(0.4885884864, 0.3919573931, Vec3(0.0, 0.0, 1.6527810236), 1, 0, 0);
    // const AtomicOrbital F_2Px(Vec3(0.0, 0.0, 1.6527810236), {F_2Px_pg1, F_2Px_pg2, F_2Px_pg3});
    //
    // const PrimitiveGaussian F_2Py_pg1(6.464803249, 0.1559162750, Vec3(0.0, 0.0, 1.6527810236), 0, 1, 0);
    // const PrimitiveGaussian F_2Py_pg2(1.502281245, 0.6076837186, Vec3(0.0, 0.0, 1.6527810236), 0, 1, 0);
    // const PrimitiveGaussian F_2Py_pg3(0.4885884864, 0.3919573931, Vec3(0.0, 0.0, 1.6527810236), 0, 1, 0);
    // const AtomicOrbital F_2Py(Vec3(0.0, 0.0, 1.6527810236), {F_2Py_pg1, F_2Py_pg2, F_2Py_pg3});
    //
    // const PrimitiveGaussian F_2Pz_pg1(6.464803249, 0.1559162750, Vec3(0.0, 0.0, 1.6527810236), 0, 0, 1);
    // const PrimitiveGaussian F_2Pz_pg2(1.502281245, 0.6076837186, Vec3(0.0, 0.0, 1.6527810236), 0, 0, 1);
    // const PrimitiveGaussian F_2Pz_pg3(0.4885884864, 0.3919573931, Vec3(0.0, 0.0, 1.6527810236), 0, 0, 1);
    // const AtomicOrbital F_2Pz(Vec3(0.0, 0.0, 1.6527810236), {F_2Pz_pg1, F_2Pz_pg2, F_2Pz_pg3});

    std::cout << "O 2Pz: " << O_2Pz.toString() << "\n";


    return 0;
}
