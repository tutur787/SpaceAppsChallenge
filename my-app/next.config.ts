import type { NextConfig } from "next";

const nextConfig: NextConfig = {
  output: 'export',
  distDir: 'embed',
  images: {
    unoptimized: true,
  },
};

export default nextConfig;
