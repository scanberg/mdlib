// The MIT License
// Copyright © 2016 Inigo Quilez
// https://www.shadertoy.com/view/4lcSRn

vec4 iCylinder( in vec3 ro, in vec3 rd, in vec3 pa, in vec3 pb, float ra ) 
{
    vec3  ba = pb - pa;
    vec3  oc = ro - pa;

    float baba = dot(ba,ba);
    float bard = dot(ba,rd);
    float baoc = dot(ba,oc);
    
    float k2 = baba            - bard*bard;
    float k1 = baba*dot(oc,rd) - baoc*bard;
    float k0 = baba*dot(oc,oc) - baoc*baoc - ra*ra*baba;
    
    float h = k1*k1 - k2*k0;
    if( h<0.0 ) return vec4(-1.0);
    h = sqrt(h);
    float t = (-k1 - h)/k2;

    float y = baoc + t*bard;
    if( y>0.0 && y<baba )
        return vec4( t, (oc+t*rd - ba*y/baba)/ra );    // side hit

    // caps
    t = (((y<0.0)?0.0:baba) - baoc)/bard;
    if( abs(k1 + k2*t) < h )
        return vec4( t, ba*sign(y)/sqrt(baba) );        // cap hit

    return vec4(-1.0);
}

// Intersect a ray with a "dashed" cylinder made of repeated capped cylinders.
vec4 iDashedCylinder(
    in vec3 ro, in vec3 rd,
    in vec3 pa, in vec3 pb,    // endpoints of the overall dashed line axis
    float ra,                  // cylinder radius
    float dashLen,             // length of one dash (along axis)
    float gapLen,              // gap length between dashes
    float dashOffset )         // phase offset (along axis)
{
    vec3  ba = pb - pa;
    float fullLen = length(ba);
    vec3  axisDir = ba / fullLen;

    float period = dashLen + gapLen;

    // project ray origin onto axis in parametric coordinate s (0..fullLen)
    float s0 = dot(ro - pa, axisDir);

    // we’ll test the dash containing s0 and its neighbors to catch intersections
    vec4 best = vec4(-1.0, 0.0, 0.0, 0.0);

    for(int i=-1; i<=1; i++){
        float sCenter = floor((s0 - dashOffset)/period)*period + float(i)*period + dashLen*0.5 + dashOffset;
        if(sCenter < 0.0 || sCenter > fullLen) continue;

        // endpoints of this dash segment in world space
        vec3 paSeg = pa + axisDir * (sCenter - dashLen*0.5);
        vec3 pbSeg = paSeg + axisDir * dashLen;

        vec4 res = iCylinder(ro, rd, paSeg, pbSeg, ra);
        if(res.x > 0.0 && (best.x < 0.0 || res.x < best.x)){
            best = res
        }
    }

    return best;
}