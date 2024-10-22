#include<iostream>
#include<string>
#include "mot_header.hh"

int main(){
    std::string str_mot_ar;
    //RNA 3D Motif Atlas sequences in forward and reverse
    str_mot_ar.append("RNA 3D Motif Atlas Hairpins forward\n");
    str_mot_ar.append(hbgsu_fw,hbgsu_fw_len);
    str_mot_ar.append("\nRNA 3D Motif Atlas Hairpins internals forward\n");
    str_mot_ar.append(ibgsu_fw,ibgsu_fw_len);
    str_mot_ar.append("\nRNA 3D Motif Atlas Bulges forward\n");
    str_mot_ar.append(bbgsu_fw,bbgsu_fw_len);
    str_mot_ar.append("\nRNA 3D Motif Atlas Hairpins reverse\n");
    str_mot_ar.append(hbgsu_rv,hbgsu_rv_len);
    str_mot_ar.append("\nRNA 3D Motif Atlas Internals reverse\n");
    str_mot_ar.append(ibgsu_rv,ibgsu_rv_len);
    str_mot_ar.append("\nRNA 3D Motif Atlas Bulges reverse\n");
    str_mot_ar.append(bbgsu_rv,bbgsu_rv_len);
    //Rfam sequences in forward and reverse
    str_mot_ar.append("\nRfam Hairpins forward\n");
    str_mot_ar.append(hrfam_fw,hrfam_fw_len);
    str_mot_ar.append("\nRfam Hairpins internals forward\n");
    str_mot_ar.append(irfam_fw,irfam_fw_len);
    str_mot_ar.append("\nRfam Bulges forward\n");
    str_mot_ar.append(brfam_fw,brfam_fw_len);
    str_mot_ar.append("\nRfam Hairpins reverse\n");
    str_mot_ar.append(hrfam_rv,hrfam_rv_len);
    str_mot_ar.append("\nRfam Internals reverse\n");
    str_mot_ar.append(irfam_rv,irfam_rv_len);
    str_mot_ar.append("\nRfam Bulges reverse\n");
    str_mot_ar.append(brfam_rv,brfam_rv_len);

    std::cout << str_mot_ar << std::endl;
}