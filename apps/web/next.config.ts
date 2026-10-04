import type { NextConfig } from "next";

const apiOrigin =
  process.env.NEXT_PUBLIC_HELIX_API_BASE_URL ||
  (process.env.NEXT_PUBLIC_VERCEL_ENV || process.env.VERCEL === "1" ? "" : "http://127.0.0.1:8000");

const cspReportOnly = [
  "default-src 'self'",
  `connect-src 'self' ${apiOrigin} http://localhost:8000 http://127.0.0.1:8000`,
  "script-src 'self' 'wasm-unsafe-eval' 'unsafe-eval'",
  "worker-src 'self' blob:",
  "style-src 'self' 'unsafe-inline'",
  "img-src 'self' data: blob:",
  "font-src 'self'",
  "media-src 'self' blob:",
  "object-src 'none'",
  "base-uri 'self'",
  "frame-ancestors 'none'",
].join("; ");

const nextConfig: NextConfig = {
  reactStrictMode: true,
  poweredByHeader: false,
  serverExternalPackages: ["plotly.js", "@mediapipe/tasks-vision"],
  async headers() {
    return [
      {
        source: "/:path*",
        headers: [
          { key: "X-Content-Type-Options", value: "nosniff" },
          { key: "Referrer-Policy", value: "no-referrer" },
          { key: "X-Frame-Options", value: "DENY" },
          { key: "Permissions-Policy", value: "camera=(self), microphone=(), geolocation=()" },
          { key: "Content-Security-Policy-Report-Only", value: cspReportOnly },
        ],
      },
    ];
  },
};

export default nextConfig;
