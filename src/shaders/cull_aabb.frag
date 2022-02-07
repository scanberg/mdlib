#version 430 core

layout(early_fragment_tests) in;
layout(std430,binding=0) buffer visiblity_buffer {
  uint visible[];
};

flat in uint id;

void main() {
    visible[id] = 1;
}